
"""
Image-based cell or nucleus segmentation using SOPA (CellPose or StarDist).

Input:
    Existing SpatialData Zarr stores.

Output:
    Reduced SpatialData Zarr stores containing segmentation results.
    => segmentation mask registered back to original Zarr store

Stages:
    1. Create reduced segmentation Zarr.
    2. Generate SOPA image patches.
    3. Run CellPose and/or StarDist segmentation per patch.
    4. Resolve segmentation boundaries (using SOPA).
    5. Register segmentation mask (shapes) to original Zarr store.

NOTE (beta):
    Do not run CellPose and StarDist simultaneously on the same
    SpatialData Zarr. Concurrent writes are not yet protected.
"""

from pathlib import Path
import shlex
import re
import subprocess
import hashlib
import json
import sys

import spatialdata as spd


# ----------------------- Global configuration ----------------------- #

configfile: "sopaSegmentImage.yaml"

RESOURCES = {"threads": 4, "mem_mb": 16000, "time": "04:00:00", "partition": "short"}
RESOURCES.update(config.get("resources", {}))

ZARR_DIR = Path(config["zarr_dir"])
SEGMENT_DIR = Path(config.get("segment_dir", "sopa.segment.dir"))

RUN_CELLPOSE = config.get("run_cellpose", False)
RUN_STARDIST = config.get("run_stardist", False)

# For beta release only!
if RUN_CELLPOSE and RUN_STARDIST:
    raise ValueError(
        "Running CellPose and StarDist together is not supported "
        "in the beta release. Enable only one method."
    )

# Shared image patching parameters
PATCH_WIDTH = config["img_patch_width"]
PATCH_OVERLAP = config["img_patch_ovlp"]

def check_owner(path, method_hash):
    actual = json.loads(Path(path).read_text()).get("hash")
    if actual != method_hash:
        raise RuntimeError(
            f"Patch cache belongs to {actual}, not {method_hash}. "
            "Rerun segmentation with the intended configuration."
        )


# ----------------------- CellPose configuration ----------------------- #

if RUN_CELLPOSE:

    CP_CONFIG = config["cp_config"]

    CP_PARAMS = {
        "patch_width": PATCH_WIDTH,
        "patch_overlap": PATCH_OVERLAP,
        "cellpose_config": CP_CONFIG,
        "cell_area": config["cp_cell_area"]
    }

    CP_HASH = hashlib.sha256(
        json.dumps(CP_PARAMS, sort_keys=True).encode()
    ).hexdigest()[:12]

    CP_MASK_KEY = config.get("cp_mask_key") or f"cellpose_{CP_HASH}"
    CP_CACHE = f"cellpose_{CP_HASH}"

    def cellpose_channels():
        channels = CP_CONFIG["channels"]

        if isinstance(channels, str):
            channels = channels.split(",")

        channels = [str(c).strip() for c in channels]

        if not channels or any(not c for c in channels):
            raise ValueError("Invalid CellPose channels.")

        return " ".join(
            f"--channels {shlex.quote(c)}"
            for c in channels
        )


# ----------------------- StarDist configuration ----------------------- #

if RUN_STARDIST:

    SD_CONFIG = config["sd_config"]
    
    if not Path(SD_CONFIG["model_path"]).is_dir():
        raise FileNotFoundError(
            f"StarDist model not found: {SD_CONFIG['model_path']}"
        )

    SD_PARAMS = {
        "patch_width": PATCH_WIDTH,
        "patch_overlap": PATCH_OVERLAP,
        "stardist_config": SD_CONFIG,
        "nucleus_area": config["sd_nucleus_area"]
    }

    SD_HASH = hashlib.sha256(
        json.dumps(SD_PARAMS, sort_keys=True).encode()
    ).hexdigest()[:12]

    SD_MASK_KEY = config.get("sd_mask_key") or f"stardist_{SD_HASH}"
    SD_CACHE = f"stardist_{SD_HASH}"


# ----------------------- Sample selection ----------------------- #

if config.get("sample_ids") is not None:
    SAMPLE_IDS = config["sample_ids"]
    
    if isinstance(SAMPLE_IDS, str):
        SAMPLE_IDS = [
            s.strip()
            for s in SAMPLE_IDS.split(",")
            if s.strip()
        ]

else:
    SAMPLE_IDS = sorted(
        p.name.removesuffix(".zarr")
        for p in ZARR_DIR.glob("*.zarr")
        if p.is_dir()
    )

if not SAMPLE_IDS:
    raise ValueError("No samples found or specified.")

missing_samples = [
    sample for sample in SAMPLE_IDS
    if not (ZARR_DIR / f"{sample}.zarr").is_dir()
]

if missing_samples:
    raise ValueError(
        f"Samples not found in {ZARR_DIR}: {', '.join(missing_samples)}"
    )

if len(SAMPLE_IDS) != len(set(SAMPLE_IDS)):
    raise ValueError("Duplicate sample IDs detected.")


def input_zarr(wildcards):
    return str(ZARR_DIR / f"{wildcards.sample}.zarr")

def segment_zarr(wildcards):
    return str(SEGMENT_DIR / f"{wildcards.sample}.zarr")


# ----------------------- Targets ----------------------- #

# Registration parameters
REGISTRATION = {}

if RUN_CELLPOSE:
    REGISTRATION["cellpose"] = {
        "source_key": "cellpose_boundaries",
        "mask_key": CP_MASK_KEY,
        "hash": CP_HASH,
        "cache": CP_CACHE,
        "parameters": CP_PARAMS
    }

if RUN_STARDIST:
    REGISTRATION["stardist"] = {
        "source_key": "stardist_boundaries",
        "mask_key": SD_MASK_KEY,
        "hash": SD_HASH,
        "cache": SD_CACHE,
        "parameters": SD_PARAMS
    }


TARGETS = []

for method, settings in REGISTRATION.items():
    TARGETS += expand(
        str(
            SEGMENT_DIR /
            "logs/{sample}_{method}Register_{mask_key}_{hash}.sentinel"
        ),
        sample=SAMPLE_IDS,
        method=method,
        mask_key=settings["mask_key"],
        hash=settings["hash"],
    )

if not TARGETS:
    TARGETS = expand(
        str(SEGMENT_DIR / "logs/{sample}_sopa_segmentSetup.sentinel"),
        sample=SAMPLE_IDS,
    )


rule full:
    input:
        TARGETS
    threads: 1
    resources:
        mem_mb=100,
        time="00:01:00",
        partition="short",
    shell:
        """
        echo "=======> All done! <======="
        """


# ----------------------- Segmentation setup ----------------------- #

rule segment_setup:
    output:
        sentinel=touch(
            str(SEGMENT_DIR / "logs/{sample}_sopa_segmentSetup.sentinel")
        )
    params:
        original_zarr=input_zarr,
        zarr=segment_zarr,
        mode="image",
        image=(
            f"--image {shlex.quote(str(config['image_key']))}"
            if config.get("image_key")
            else ""
        ),
    resources:
        **RESOURCES,
    log:
        str(SEGMENT_DIR / "logs/{sample}_sopa_segmentSetup.log")
    shell:
        """
        python -m spatialhub.scripts.sopaSegmentSetup \
            --input-zarr {params.original_zarr:q} \
            --output-zarr {params.zarr:q} \
            --segmentation-mode {params.mode} \
            {params.image} \
            --overwrite \
            > {log:q} 2>&1

        test -d {params.zarr:q}
        """


# ----------------------- Image patching ----------------------- #

checkpoint sopa_image_patch:
    input:
        str(SEGMENT_DIR / "logs/{sample}_sopa_segmentSetup.sentinel")
    output:
        patch_count=str(SEGMENT_DIR / "patches/{sample}_image.count")
    params:
        zarr=segment_zarr,
        width=PATCH_WIDTH,
        overlap=PATCH_OVERLAP,
    resources:
        **RESOURCES,
    log:
        str(SEGMENT_DIR / "logs/{sample}_patchImage.log")
    run:
        command = [
            "sopa", "patchify", "image", params.zarr,
            "--patch-width-pixel", str(params.width),
            "--patch-overlap-pixel", str(params.overlap),
        ]

        with open(log[0], "w") as handle:
            subprocess.run(
                command,
                stdout=handle,
                stderr=subprocess.STDOUT,
                check=True,
            )

        with open(log[0]) as handle:
            content = handle.read()

        matches = re.findall(
            r"Added (\d+) patch\(es\) to sdata\['image_patches'\]",
            content
        )

        if not matches:
            raise ValueError(
                f"Could not determine patch count from {log[0]}"
            )

        n_patches = int(matches[-1])

        if n_patches < 1:
            raise ValueError("No image patches were generated.")

        with open(output.patch_count, "w") as handle:
            handle.write(f"{n_patches}\n")



# ----------------------- Patch output parsing ----------------------- #

def image_settings(wc):
    settings = REGISTRATION[wc.method]

    for key in ("hash", "mask_key"):
        if getattr(wc, key, settings[key]) != settings[key]:
            raise ValueError(f"Invalid {key} for {wc.method}")

    return settings

def get_patch_count(wildcards):
    ckpt = checkpoints.sopa_image_patch.get(
        sample=wildcards.sample
    )

    with open(ckpt.output.patch_count) as handle:
        n_patches = int(handle.read().strip())

    if n_patches < 1:
        raise ValueError(
            f"No patches found for {wildcards.sample}"
        )

    return n_patches

def image_patch_outputs(wc):
    settings = image_settings(wc)

    return expand(
        str(
            SEGMENT_DIR /
            f"{wc.sample}.zarr/.sopa_cache/{settings['cache']}/{{patch}}.parquet"
        ),
        patch=range(get_patch_count(wc)),
    )


# Track the active configuration to invalidate stale resolution and registration
rule track_image:
    input:
        str(SEGMENT_DIR / "patches/{sample}_image.count")
    output:
        current=str(
            SEGMENT_DIR / "patches/{sample}_{method}_current.json"
        )
    params:
        hash=lambda wc: image_settings(wc)["hash"]
    wildcard_constraints:
        method="cellpose|stardist"
    resources:
        **RESOURCES
    run:
        Path(output.current).write_text(
            json.dumps({"hash": params.hash}) + "\n"
        )


# ----------------------- CellPose segmentation ----------------------- #

if RUN_CELLPOSE:

    rule sopa_segment_cellpose:
        input:
            str(SEGMENT_DIR / "patches/{sample}_image.count")
        output:
            parquet=str(
                SEGMENT_DIR /
                f"{{sample}}.zarr/.sopa_cache/{CP_CACHE}/{{patch}}.parquet"
            )
        params:
            zarr=segment_zarr,
            channels=cellpose_channels(),
            model=CP_CONFIG["model"],
            diameter=CP_CONFIG["cell_diameter"],
            min_area=config["cp_cell_area"],
            flow_threshold=CP_CONFIG["flow_thresh"],
            prob_threshold=CP_CONFIG["prob_thresh"],
            clip_limit=CP_CONFIG["clip_limit"],
            gaussian_sigma=CP_CONFIG["gaussian_sigma"],
            clahe_kernel=(
                f"--clahe-kernel-size {CP_CONFIG['clahe_kernel_size']}"
                if CP_CONFIG.get("clahe_kernel_size") is not None
                and str(CP_CONFIG["clahe_kernel_size"]).strip()
                else ""
            ),
            cache=CP_CACHE,
        wildcard_constraints:
            patch=r"\d+"
        resources:
            **RESOURCES
        log:
            str(
                SEGMENT_DIR /
                f"logs/{{sample}}_sopaSegment_{CP_HASH}_{{patch}}.log"
            )
        shell:
            """
            sopa segmentation cellpose {params.zarr:q} \
                {params.channels} \
                --pretrained-model {params.model:q} \
                --diameter {params.diameter} \
                --min-area {params.min_area} \
                --flow-threshold {params.flow_threshold} \
                --cellprob-threshold {params.prob_threshold} \
                --clip-limit {params.clip_limit} \
                --gaussian-sigma {params.gaussian_sigma} \
                {params.clahe_kernel} \
                --patch-index {wildcards.patch} \
                --cache-dir-name {params.cache:q} \
                > {log:q} 2>&1
            """


# ----------------------- StarDist segmentation ----------------------- #

if RUN_STARDIST:
    
    rule sopa_segment_stardist:
        input:
            str(SEGMENT_DIR / "patches/{sample}_image.count")
        output:
            parquet=str(
                SEGMENT_DIR /
                f"{{sample}}.zarr/.sopa_cache/{SD_CACHE}/{{patch}}.parquet"
            )
        params:
            zarr=segment_zarr,
            channels=SD_CONFIG["channels"],
            model=SD_CONFIG["model_path"],
            prob_threshold=SD_CONFIG["prob_thresh"],
            nms_threshold=SD_CONFIG["nms_thresh"],
            min_area=config["sd_nucleus_area"],
            clip_limit=SD_CONFIG["clip_limit"],
            gaussian_sigma=SD_CONFIG["gaussian_sigma"],
            clahe_kernel=(
                f"--clahe-kernel-size {SD_CONFIG['clahe_kernel_size']}"
                if SD_CONFIG.get("clahe_kernel_size") is not None
                and str(SD_CONFIG["clahe_kernel_size"]).strip()
                else ""
            ),
            cache=SD_CACHE,
        wildcard_constraints:
            patch=r"\d+"
        resources:
            **RESOURCES
        log:
            str(
                SEGMENT_DIR /
                f"logs/{{sample}}_sopaSegment_{SD_HASH}_{{patch}}.log"
            )
        shell:
            """
            sopa segmentation stardist {params.zarr:q} \
                --channels {params.channels:q} \
                --local-model {params.model:q} \
                --prob-thresh {params.prob_threshold} \
                --nms-thresh {params.nms_threshold} \
                --min-area {params.min_area} \
                --clip-limit {params.clip_limit} \
                --gaussian-sigma {params.gaussian_sigma} \
                {params.clahe_kernel} \
                --patch-index {wildcards.patch} \
                --cache-dir-name {params.cache:q} \
                > {log:q} 2>&1
            """


# ----------------------- Image-based mask resolution ----------------------- #

rule image_completion:
    input:
        patches=image_patch_outputs,
        current=str(
            SEGMENT_DIR / "patches/{sample}_{method}_current.json"
        )
    output:
        completed=str(
            SEGMENT_DIR /
            "configs/{sample}/{method}Complete_{hash}.json"
        )
    params:
        settings=image_settings
    wildcard_constraints:
        method="cellpose|stardist",
        hash=r"[0-9a-f]{12}"
    resources:
        **RESOURCES
    run:
        settings = params.settings
        check_owner(input.current, settings["hash"])

        state = {
            "hash": settings["hash"],
            "status": "complete",
            "parameters": settings["parameters"],
        }

        Path(output.completed).write_text(
            json.dumps(state, indent=2, sort_keys=True) + "\n"
        )


rule sopa_resolve_image:
    input:
        completed=str(
            SEGMENT_DIR /
            "configs/{sample}/{method}Complete_{hash}.json"
        ),
        current=str(
            SEGMENT_DIR / "patches/{sample}_{method}_current.json"
        )
    output:
        sentinel=touch(
            str(
                SEGMENT_DIR /
                "patches/{sample}/{method}Resolve_{hash}.sentinel"
            )
        )
    params:
        zarr=segment_zarr,
        settings=image_settings
    wildcard_constraints:
        method="cellpose|stardist",
        hash=r"[0-9a-f]{12}"
    resources:
        **RESOURCES
    log:
        str(
            SEGMENT_DIR /
            "logs/{sample}_{method}Resolve_{hash}.log"
        )
    run:
        settings = params.settings
        check_owner(input.current, settings["hash"])

        completion = json.loads(Path(input.completed).read_text())

        if completion != {
            "hash": settings["hash"],
            "status": "complete",
            "parameters": settings["parameters"],
        }:
            raise ValueError("Image segmentation completion record is invalid.")

        command = [
            "sopa", "resolve", wildcards.method, params.zarr,
            "--cache-dir-name", settings["cache"],
        ]

        with open(log[0], "w") as handle:
            subprocess.run(
                command,
                stdout=handle,
                stderr=subprocess.STDOUT,
                check=True,
            )

        sdata = spd.read_zarr(params.zarr)

        if settings["source_key"] not in sdata.shapes:
            raise RuntimeError(
                f"Resolution did not produce '{settings['source_key']}'."
            )


# ----------------------- Registration to source zarr ----------------------- #

def registration_source(wc):
    settings = image_settings(wc)

    return str(
        SEGMENT_DIR /
        f"patches/{wc.sample}/{wc.method}Resolve_"
        f"{settings['hash']}.sentinel"
    )

rule sopa_register_mask:
    input:
        sentinel=registration_source,
        current=str(
            SEGMENT_DIR / "patches/{sample}_{method}_current.json"
        )
    output:
        sentinel=touch(
            str(
                SEGMENT_DIR /
                "logs/{sample}_{method}Register_{mask_key}_{hash}.sentinel"
            )
        )
    params:
        zarr=segment_zarr,
        original_zarr=input_zarr,
        source_key=lambda wc: REGISTRATION[wc.method]["source_key"],
        mask_key=lambda wc: REGISTRATION[wc.method]["mask_key"],
        image_key=config.get("image_key", "image"),
        overwrite=config.get("overwrite_mask", False),
    wildcard_constraints:
        method="cellpose|stardist",
        hash=r"[0-9a-f]{12}"
    resources:
        **RESOURCES
    log:
        str(
            SEGMENT_DIR /
            "logs/{sample}_{method}Register_{mask_key}_{hash}.log"
        )
    run:
        check_owner(
            input.current,
            REGISTRATION[wildcards.method]["hash"],
        )
        
        command = [
            sys.executable,
            "-m", "spatialhub.scripts.sopaRegisterMask",
            "--input-zarr", params.zarr,
            "--output-zarr", params.original_zarr,
            "--source-key", params.source_key,
            "--mask-key", params.mask_key,
            "--image-key", params.image_key,
        ]

        if params.overwrite:
            command.append("--overwrite")

        with open(log[0], "w") as handle:
            subprocess.run(
                command,
                stdout=handle,
                stderr=subprocess.STDOUT,
                check=True,
            )