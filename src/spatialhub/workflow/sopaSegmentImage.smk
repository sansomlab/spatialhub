
"""
Image-based cell segmentation using SOPA/Cellpose.

Input:
    Existing SpatialData Zarr stores.

Output:
    Reduced SpatialData Zarr stores containing segmentation results.

Stages:
    1. Create reduced segmentation Zarr.
    2. Generate SOPA image patches.
    3. Run Cellpose segmentation.
    4. Register segmentation boundaries.
"""

from pathlib import Path
import shlex
import re
import subprocess
import hashlib
import json


# ----------------------- Configuration ----------------------- #

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

CP_MODEL = config["cp_model"]
CP_MASK_KEY = config["cp_mask_key"]
RUN_CELLPOSE = config.get("run_cellpose", True)

CP_PARAMS = {
    "patch_width": config["img_patch_width"],
    "patch_overlap": config["img_patch_ovlp"],
    "model": config["cp_model"],
    "mask_key": config["cp_mask_key"],
    "channels": config["cp_channels"],
    "diameter": config["cp_cell_diameter"],
    "min_area": config["cp_cell_area"],
    "flow_threshold": config["cp_flow_thresh"],
    "prob_threshold": config["cp_prob_thresh"],
}

CP_HASH = hashlib.sha256(
    json.dumps(CP_PARAMS, sort_keys=True).encode()
).hexdigest()[:12]

CP_CACHE = f"cellpose_{CP_HASH}"


# ----------------------- Sample selection ----------------------- #

if config.get("sample_ids") is not None:
    SAMPLE_IDS = config["sample_ids"]
else:
    SAMPLE_IDS = sorted(
        p.name.removesuffix(".zarr")
        for p in ZARR_DIR.glob("*.zarr")
        if p.is_dir()
    )

if not SAMPLE_IDS:
    raise ValueError("No samples found or specified.")

if len(SAMPLE_IDS) != len(set(SAMPLE_IDS)):
    raise ValueError("Duplicate sample IDs detected.")


def input_zarr(wildcards):
    return str(ZARR_DIR / f"{wildcards.sample}.zarr")


def segment_zarr(wildcards):
    return str(SEGMENT_DIR / f"{wildcards.sample}.zarr")


def cellpose_channels():
    channels = config["cp_channels"]

    if isinstance(channels, str):
        channels = channels.split(",")

    return " ".join(
        f"--channels {shlex.quote(str(channel).strip())}"
        for channel in channels
    )


# ----------------------- Targets ----------------------- #

rule full:
    input:
        expand(
            str(
                SEGMENT_DIR /
                f"patches/{{sample}}/cellposeResolve_{CP_MASK_KEY}_{CP_HASH}.sentinel"
            ),
            sample=SAMPLE_IDS,
        ) if RUN_CELLPOSE else expand(
            str(SEGMENT_DIR / "logs/{sample}_sopa_segmentSetup.sentinel"),
            sample=SAMPLE_IDS,
        )
    params:
        lock=LOCK,
    resources:
        **RESOURCES,
    shell:
        """
        echo "=======> All done! <======="
        """


# ----------------------- Segmentation setup ----------------------- #

#RESOURCES["mem_mb"] = config["sdata_mem"] 

rule segment_setup:
    input:
        zarr=input_zarr
    output:
        sentinel=touch(
            str(SEGMENT_DIR / "logs/{sample}_sopa_segmentSetup.sentinel")
        )
    params:
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
            --input-zarr {input.zarr:q} \
            --output-zarr {params.zarr:q} \
            --segmentation-mode {params.mode} \
            {params.image} \
            {params.points} \
            {params.prior} \
            --overwrite \
            > {log:q} 2>&1
        """


# ----------------------- Cellpose patching ----------------------- #

#RESOURCES["mem_mb"] = config["cp_mem"]
#RESOURCES["threads"] = config["cp_cpu"] 

checkpoint sopa_image_patch:
    input:
        str(SEGMENT_DIR / "logs/{sample}_sopa_segmentSetup.sentinel")
    output:
        sentinel=touch(
            str(SEGMENT_DIR / "logs/{sample}_sopa_imagePatch.sentinel")
        ),
        patch_count=str(SEGMENT_DIR / "patches/{sample}_image.count")
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

        with open(output.patch_count, "w") as handle:
            handle.write(f"{n_patches}\n")



def image_patch_outputs(wildcards):

    ckpt = checkpoints.sopa_image_patch.get(
        sample=wildcards.sample
    )

    with open(ckpt.output.patch_count) as handle:
        n_patches = int(handle.read().strip())

    if n_patches < 1:
        raise ValueError(
            f"No patches found for {wildcards.sample}"
        )

    return expand(
        str(
            SEGMENT_DIR /
            f"patches/{{sample}}/cellposeSegment_{CP_MASK_KEY}_{CP_HASH}_{{patch}}.sentinel"
        ),
        sample=wildcards.sample,
        patch=range(n_patches),
    )


# ----------------------- Cellpose segmentation ----------------------- #

rule sopa_segment_cellpose:
    input:
        str(SEGMENT_DIR / "logs/{sample}_sopa_imagePatch.sentinel")
    output:
        sentinel=touch(
            str(SEGMENT_DIR / f"patches/{{sample}}/cellposeSegment_{CP_MASK_KEY}_{CP_HASH}_{{patch}}.sentinel")
        )
    params:
        zarr=segment_zarr,
        channels=cellpose_channels(),
        model=CP_MODEL,
        diameter=config["cp_cell_diameter"],
        min_area=config["cp_cell_area"],
        flow_threshold=config["cp_flow_thresh"],
        prob_threshold=config["cp_prob_thresh"],
        cache=CP_CACHE,
    wildcard_constraints:
        patch=r"\d+"
    resources:
        **RESOURCES,
    log:
        str(
            SEGMENT_DIR /
            f"logs/{{sample}}_sopaSegment_{CP_MASK_KEY}_{CP_HASH}_{{patch}}.log"
        )
    shell:
        """
        sopa segmentation cellpose {params.zarr} \
            {params.channels} \
            --pretrained-model {params.model:q} \
            --diameter {params.diameter} \
            --min-area {params.min_area} \
            --flow-threshold {params.flow_threshold} \
            --cellprob-threshold {params.prob_threshold} \
            --patch-index {wildcards.patch} \
            --cache-dir-name {params.cache:q} \
            > {log:q} 2>&1
        """


rule sopa_resolve_cellpose:
    input:
        image_patch_outputs
    output:
        sentinel=touch(
            str(
                SEGMENT_DIR /
                f"patches/{{sample}}/cellposeResolve_{CP_MASK_KEY}_{CP_HASH}.sentinel"
            )
        )
    params:
        zarr=segment_zarr,
        cache=CP_CACHE,
        mask_key=CP_MASK_KEY,
    resources:
        **RESOURCES,
    log:
        str(
            SEGMENT_DIR /
            f"logs/{{sample}}_sopaResolve_{CP_MASK_KEY}_{CP_HASH}.log"
        )
    run:
        zarr = Path(params.zarr)
        source = zarr / "shapes" / "cellpose_boundaries"
        destination = zarr / "shapes" / params.mask_key

        if source != destination and destination.exists():
            raise FileExistsError(
                f"Segmentation already exists: {destination}"
            )

        command = [
            "sopa", "resolve", "cellpose", str(zarr),
            "--cache-dir-name", params.cache,
        ]

        with open(log[0], "w") as handle:
            subprocess.run(
                command,
                stdout=handle,
                stderr=subprocess.STDOUT,
                check=True,
            )

        if not source.exists():
            raise FileNotFoundError(
                f"Resolved segmentation not found: {source}"
            )

        if source != destination:
            source.rename(destination)

        with open(log[0], "a") as handle:
            handle.write(
                f"\nRegistered segmentation: {destination}\n"
            )

