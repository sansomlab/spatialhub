
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

# [NOTE] this key is set by the CLI, not the config file!
if "lock" in config:
    LOCK = config["lock"].lower()
else:
    LOCK = "false"

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

        return " ".join(
            f"--channels {shlex.quote(str(c).strip())}"
            for c in channels
        )


# ----------------------- Stardist configuration ----------------------- #

if RUN_STARDIST:

    SD_CONFIG = config["sd_config"]
    SD_MASK_KEY = config["sd_mask_key"]

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

TARGETS = []

if RUN_CELLPOSE:
    TARGETS += expand(
        str(
            SEGMENT_DIR /
            f"patches/{{sample}}/cellposeResolve_{CP_MASK_KEY}_{CP_HASH}.sentinel"
        ),
        sample=SAMPLE_IDS,
    )

if RUN_STARDIST:
    TARGETS += expand(
        str(
            SEGMENT_DIR /
            f"patches/{{sample}}/stardistResolve_{SD_MASK_KEY}_{SD_HASH}.sentinel"
        ),
        sample=SAMPLE_IDS,
    )

if not TARGETS:
    TARGETS = expand(
        str(SEGMENT_DIR / "logs/{sample}_sopa_segmentSetup.sentinel"),
        sample=SAMPLE_IDS,
    )


rule full:
    input:
        TARGETS
    params:
        lock=LOCK
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
        points=(
            f"--points {shlex.quote(str(config['points_key']))}"
            if config.get("points_key")
            else ""
        ),
        prior=(
            f"--prior {shlex.quote(str(config['prior_segment_key']))}"
            if config.get("prior_segment_key")
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
            {params.points} \
            {params.prior} \
            --overwrite \
            > {log:q} 2>&1
        """


# ----------------------- Cellpose patching ----------------------- #

checkpoint sopa_image_patch:
    input:
        str(SEGMENT_DIR / "logs/{sample}_sopa_segmentSetup.sentinel")
    output:
        patch_indices=str(SEGMENT_DIR / "patches/{sample}_image.indices")
    params:
        zarr=segment_zarr,
        width=config["img_patch_width"],
        overlap=config["img_patch_ovlp"],
    resources:
        **RESOURCES,
    log:
        str(SEGMENT_DIR / "logs/{sample}_patchIamge.log")
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

        with open(output.patch_indices, "w") as handle:
            handle.write(f"{n_patches}\n")



# ----------------------- Patch output parsing ----------------------- #

def get_patch_indices(wildcards):
    ckpt = checkpoints.sopa_image_patch.get(
        sample=wildcards.sample
    )

    with open(ckpt.output.patch_indices) as handle:
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
                f"patches/{{sample}}/cellposeSegment_{CP_MASK_KEY}_{CP_HASH}_{{patch}}.sentinel"
            ),
            sample=wildcards.sample,
            patch=range(get_patch_indices(wildcards)),
        )


if RUN_STARDIST:

    def stardist_patch_outputs(wildcards):
        return expand(
            str(
                SEGMENT_DIR /
                f"patches/{{sample}}/stardistSegment_{SD_MASK_KEY}_{SD_HASH}_{{patch}}.sentinel"
            ),
            sample=wildcards.sample,
            patch=range(get_patch_indices(wildcards)),
        )


# ----------------------- Cellpose segmentation ----------------------- #

if RUN_CELLPOSE:

    rule sopa_segment_cellpose:
        input:
            str(SEGMENT_DIR / "patches/{sample}_image.indices")
        output:
            sentinel=touch(
                str(
                    SEGMENT_DIR /
                    f"patches/{{sample}}/cellposeSegment_{CP_MASK_KEY}_{CP_HASH}_{{patch}}.sentinel"
                )
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
            sentinel=touch(
                str(
                    SEGMENT_DIR /
                    f"patches/{{sample}}/stardistSegment_{SD_MASK_KEY}_{SD_HASH}_{{patch}}.sentinel"
                )
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
            )
        params:
            zarr=segment_zarr,
            original_zarr=input_zarr,
            cache=CP_CACHE,
            mask_key=CP_MASK_KEY,
            image_key=config.get("image_key", "image"),
        resources:
            **RESOURCES
        log:
            str(
                SEGMENT_DIR /
                f"logs/{{sample}}_sopaResolve_{CP_MASK_KEY}_{CP_HASH}.log"
            )
        run:
            import sys

            command = [
                "sopa", "resolve", "cellpose", params.zarr,
                "--cache-dir-name", params.cache,
            ]

            with open(log[0], "w") as handle:
                subprocess.run(
                    command,
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                    check=True,
                )

            command = [
                sys.executable,
                "-m", "spatialhub.scripts.sopaRegisterMask",
                "--input-zarr", params.zarr,
                "--output-zarr", params.original_zarr,
                "--source-key", "cellpose_boundaries",
                "--mask-key", params.mask_key,
                "--image-key", params.image_key,
            ]

            with open(log[0], "a") as handle:
                subprocess.run(
                    command,
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                    check=True,
                )


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
            )
        params:
            zarr=segment_zarr,
            original_zarr=input_zarr,
            cache=SD_CACHE,
            mask_key=SD_MASK_KEY,
            image_key=config.get("image_key", "image"),
        resources:
            **RESOURCES
        log:
            str(
                SEGMENT_DIR /
                f"logs/{{sample}}_sopaResolve_{SD_MASK_KEY}_{SD_HASH}.log"
            )
        run:
            import sys

            command = [
                "sopa", "resolve", "stardist", params.zarr,
                "--cache-dir-name", params.cache,
            ]

            with open(log[0], "w") as handle:
                subprocess.run(
                    command,
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                    check=True,
                )

            command = [
                sys.executable,
                "-m", "spatialhub.scripts.sopaRegisterMask",
                "--input-zarr", params.zarr,
                "--output-zarr", params.original_zarr,
                "--source-key", "stardist_boundaries",
                "--mask-key", params.mask_key,
                "--image-key", params.image_key,
            ]

            with open(log[0], "a") as handle:
                subprocess.run(
                    command,
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                    check=True,
                )
