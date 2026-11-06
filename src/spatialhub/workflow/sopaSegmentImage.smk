
"""
Image-based cell or nucleus segmentation using SOPA (Cellpose or StarDist).

Input:
    Existing SpatialData Zarr stores.

Output:
    Reduced SpatialData Zarr stores containing segmentation results.

Stages:
    1. Create reduced segmentation Zarr.
    2. Generate SOPA image patches.
    3. Run Cellpose and/or StarDist segmentation per patch.
    4. Resolve segmentation boundaries.

NOTE (beta):
    Do not run Cellpose and StarDist simultaneously on the same
    SpatialData Zarr. Concurrent writes are not yet protected.
"""

from pathlib import Path
import shlex
import re
import subprocess
import hashlib
import json


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
        "Running Cellpose and StarDist together is not supported "
        "in the beta release. Enable only one method."
    )

# Shared image patching parameters
PATCH_WIDTH = config["img_patch_width"]
PATCH_OVERLAP = config["img_patch_ovlp"]


# ----------------------- Cellpose configuration ----------------------- #

if RUN_CELLPOSE:

    CP_CONFIG = config["cp_config"]
    CP_MASK_KEY = config["cp_mask_key"]

    CP_PARAMS = {
        "patch_width": PATCH_WIDTH,
        "patch_overlap": PATCH_OVERLAP,
        "cellpose_config": CP_CONFIG,
        "cell_area": config["cp_cell_area"],
        "mask_key": CP_MASK_KEY,
    }

    CP_HASH = hashlib.sha256(
        json.dumps(CP_PARAMS, sort_keys=True).encode()
    ).hexdigest()[:12]

    CP_CACHE = f"cellpose_{CP_HASH}"

    def cellpose_channels():
        channels = CP_CONFIG["channels"]

        if isinstance(channels, str):
            channels = channels.split(",")

        channels = [str(c).strip() for c in channels]

        if not channels or any(not c for c in channels):
            raise ValueError("Invalid Cellpose channels.")

        return " ".join(
            f"--channels {shlex.quote(c)}"
            for c in channels
        )


# ----------------------- Stardist configuration ----------------------- #

if RUN_STARDIST:

    SD_CONFIG = config["sd_config"]
    SD_MASK_KEY = config["sd_mask_key"]

    if not Path(SD_CONFIG["model_path"]).is_dir():
        raise FileNotFoundError(
            f"StarDist model not found: {SD_CONFIG['model_path']}"
        )

    SD_PARAMS = {
        "patch_width": PATCH_WIDTH,
        "patch_overlap": PATCH_OVERLAP,
        "stardist_config": SD_CONFIG,
        "nucleus_area": config["sd_nucleus_area"],
        "mask_key": SD_MASK_KEY,
    }

    SD_HASH = hashlib.sha256(
        json.dumps(SD_PARAMS, sort_keys=True).encode()
    ).hexdigest()[:12]

    SD_CACHE = f"stardist_{SD_HASH}"


# ----------------------- Sample selection ----------------------- #

if config.get("sample_ids") is not None:
    SAMPLE_IDS = config["sample_ids"]
    
    if isinstance(SAMPLE_IDS, str):
        SAMPLE_IDS = [s.strip() for s in SAMPLE_IDS.split(",")]

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
    }

if RUN_STARDIST:
    REGISTRATION["stardist"] = {
        "source_key": "stardist_boundaries",
        "mask_key": SD_MASK_KEY,
        "hash": SD_HASH,
    }


TARGETS = []

for method, settings in REGISTRATION.items():
    TARGETS += expand(
        str(
            SEGMENT_DIR /
            "logs/{sample}_sopaRegister_{method}_{hash}.sentinel"
        ),
        sample=SAMPLE_IDS,
        method=method,
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
    resources:
        **RESOURCES
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
        mode=config.get("segmentation_mode", "image"),
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


# ----------------------- Cellpose patching ----------------------- #

checkpoint sopa_image_patch:
    input:
        str(SEGMENT_DIR / "logs/{sample}_sopa_segmentSetup.sentinel")
    output:
        patch_count=str(SEGMENT_DIR / "patches/{sample}_image.count")
    params:
        zarr=segment_zarr,
        width=config["img_patch_width"],
        overlap=config["img_patch_ovlp"],
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


if RUN_CELLPOSE:
    def cellpose_patch_outputs(wildcards):
        return expand(
            str(
                SEGMENT_DIR /
                f"{{sample}}.zarr/.sopa_cache/{CP_CACHE}/{{patch}}.parquet"
            ),
            sample=wildcards.sample,
            patch=range(get_patch_count(wildcards)),
        )


if RUN_STARDIST:
    def stardist_patch_outputs(wildcards):
        return expand(
            str(
                SEGMENT_DIR /
                f"{{sample}}.zarr/.sopa_cache/{SD_CACHE}/{{patch}}.parquet"
            ),
            sample=wildcards.sample,
            patch=range(get_patch_count(wildcards)),
        )


# ----------------------- Cellpose segmentation ----------------------- #

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
                f"logs/{{sample}}_sopaSegment_{CP_MASK_KEY}_{CP_HASH}_{{patch}}.log"
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
                f"logs/{{sample}}_sopaSegment_{SD_MASK_KEY}_{SD_HASH}_{{patch}}.log"
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


# ----------------------- Cellpose resolution ----------------------- #

if RUN_CELLPOSE:

    rule sopa_resolve_cellpose:
        input:
            cellpose_patch_outputs
        output:
            sentinel=touch(
                str(
                    SEGMENT_DIR /
                    f"patches/{{sample}}/cellposeResolve_{CP_MASK_KEY}_{CP_HASH}.sentinel"
                )
            ),
        params:
            zarr=segment_zarr,
            cache=CP_CACHE,
        resources:
            **RESOURCES
        log:
            str(
                SEGMENT_DIR /
                f"logs/{{sample}}_sopaResolve_{CP_MASK_KEY}_{CP_HASH}.log"
            )
        shell:
            """
            sopa resolve cellpose {params.zarr:q} \
                --cache-dir-name {params.cache:q} \
                > {log:q} 2>&1
            """


# ----------------------- StarDist resolution ----------------------- #

if RUN_STARDIST:

    rule sopa_resolve_stardist:
        input:
            stardist_patch_outputs
        output:
            sentinel=touch(
                str(
                    SEGMENT_DIR /
                    f"patches/{{sample}}/stardistResolve_{SD_MASK_KEY}_{SD_HASH}.sentinel"
                )
            ),
        params:
            zarr=segment_zarr,
            cache=SD_CACHE,
        resources:
            **RESOURCES
        log:
            str(
                SEGMENT_DIR /
                f"logs/{{sample}}_sopaResolve_{SD_MASK_KEY}_{SD_HASH}.log"
            )
        shell:
            """
            sopa resolve stardist {params.zarr:q} \
                --cache-dir-name {params.cache:q} \
                > {log:q} 2>&1
            """


# ----------------------- Registration to source zarr ----------------------- #

def registration_source(wildcards):
    method = wildcards.method

    if method not in REGISTRATION:
        raise ValueError(f"Unknown segmentation method: {method}")

    expected_hash = REGISTRATION[method]["hash"]

    if wildcards.hash != expected_hash:
        raise ValueError(
            f"Configuration hash mismatch for {method}: "
            f"{wildcards.hash} != {expected_hash}"
        )

    return str(
        SEGMENT_DIR /
        f"patches/{wildcards.sample}/"
        f"{method}Resolve_{REGISTRATION[method]['mask_key']}_{expected_hash}.sentinel"
    )


rule sopa_register_mask:
    input:
        sentinel=registration_source
    output:
        sentinel=touch(
            str(
                SEGMENT_DIR /
                "logs/{sample}_sopaRegister_{method}_{hash}.sentinel"
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
            "logs/{sample}_sopaRegister_{method}_{hash}.log"
        )
    run:
        import sys

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