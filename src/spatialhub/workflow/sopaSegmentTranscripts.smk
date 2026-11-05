
"""
Transcripts-based cell segmentation using SOPA (baysor/proseg).

Input:
    Existing SpatialData Zarr stores.

Output:
    Reduced SpatialData Zarr stores containing segmentation results.

Stages:
    1. Create reduced segmentation Zarr.
    2. Generate SOPA transcripts patches.
    3. Run baysor and/or proseg segmentation.
    4. Register segmentation boundaries.
"""

from pathlib import Path
import shlex
import re
import subprocess
import hashlib
import json
import spatialdata as spd


# ----------------------- Configuration ----------------------- #

configfile: "sopaSegmentTranscripts.yaml"

# [NOTE] this key is set by the CLI, not the config file!
if "lock" in config:
    LOCK = config["lock"].lower()
else:
    LOCK = "false"

RESOURCES = {"threads": 4, "mem_mb": 16000, "time": "04:00:00", "partition": "short"}
RESOURCES.update(config.get("resources", {}))

ZARR_DIR = Path(config["zarr_dir"])
SEGMENT_DIR = Path(config.get("segment_dir", "sopa.segment.dir"))

RUN_BAYSOR = config.get("run_baysor", True)
RUN_PROSEG = config.get("run_proseg", False)

BS_PARAMS = {
    "patch_width": config["tx_patch_width"],
    "patch_overlap": config["tx_patch_ovlp"],
    "baysor_config": config["bs_config"],
    "baysor_area": config["bs_cell_area"],
    "mask_key": config["bs_mask_key"],
}

BS_HASH = hashlib.sha256(
    json.dumps(BS_PARAMS, sort_keys=True).encode()
).hexdigest()[:12]

BS_CACHE = f"baysor_{BS_HASH}"
BS_MASK_KEY = config["bs_mask_key"]


# ----------------------- Sample selection ----------------------- #

if config.get("sample_ids") is not None:
    SAMPLE_IDS = config["sample_ids"]
    # add check that this is a valid ID in the zarr path, with informative error
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

rule full:
    input:
        expand(
            str(
                SEGMENT_DIR /
                f"patches/{{sample}}/baysorResolve_{BS_MASK_KEY}_{BS_HASH}.sentinel"
            ),
            sample=SAMPLE_IDS,
        ) if RUN_BAYSOR else expand(
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
        mode=config.get("segmentation_mode", "transcripts"),
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


# ----------------------- Transcripts patching ----------------------- #

#RESOURCES["mem_mb"] = config["bs_mem"]
#RESOURCES["threads"] = config["bs_cpu"] 

checkpoint sopa_transcripts_patch:
    input:
        str(SEGMENT_DIR / "logs/{sample}_sopa_segmentSetup.sentinel")
    output:
        sentinel=touch(
            str(SEGMENT_DIR / "logs/{sample}_sopa_transcriptsPatch.sentinel")
        ),
        patch_indices=str(
            SEGMENT_DIR / "patches/{sample}_transcripts.indices"
        )
    params:
        zarr=segment_zarr,
        width=config["tx_patch_width"],
        overlap=config["tx_patch_ovlp"],
        prior=config["prior_segment_key"]
    resources:
        **RESOURCES,
    log:
        str(SEGMENT_DIR / "logs/{sample}_patchTranscripts.log")
    run:
        command = [
            "sopa", "patchify", "transcripts", params.zarr,
            "--patch-width-microns", str(params.width),
            "--patch-overlap-microns", str(params.overlap),
            "--prior-shapes-key", str(params.prior),
        ]

        with open(log[0], "w") as handle:
            subprocess.run(
                command,
                stdout=handle,
                stderr=subprocess.STDOUT,
                check=True,
            )

        # Retrieve actual patch indices (which may not be consecutive)
        sdata = spd.read_zarr(params.zarr)

        patch_indices = sdata["transcripts_patches"].index.tolist()

        if not patch_indices:
            raise ValueError("No transcript patches were generated.")

        with open(output.patch_indices, "w") as handle:
            handle.write(",".join(map(str, patch_indices)))



def transcripts_patch_outputs(wildcards):

    ckpt = checkpoints.sopa_transcripts_patch.get(
        sample=wildcards.sample
    )

    with open(ckpt.output.patch_indices) as handle:
        patch_indices = [
            int(p) for p in handle.read().strip().split(",")
        ]

    return expand(
        str(
            SEGMENT_DIR /
            f"patches/{{sample}}/baysorSegment_{BS_MASK_KEY}_{BS_HASH}_{{patch}}.sentinel"
        ),
        sample=wildcards.sample,
        patch=patch_indices,
    )


# ----------------------- Baysor segmentation ----------------------- #

rule sopa_segment_baysor:
    input:
        str(SEGMENT_DIR / "logs/{sample}_sopa_transcriptsPatch.sentinel")
    output:
        sentinel=touch(
            str(SEGMENT_DIR / f"patches/{{sample}}/baysorSegment_{BS_MASK_KEY}_{BS_HASH}_{{patch}}.sentinel")
        )
    params:
        zarr=segment_zarr,
        config=repr(config["bs_config"]),
        cache=BS_CACHE,
        min_area=config["bs_cell_area"]
    wildcard_constraints:
        patch=r"\d+"
    resources:
        **RESOURCES,
    log:
        str(
            SEGMENT_DIR /
            f"logs/{{sample}}_sopaSegment_{BS_MASK_KEY}_{BS_HASH}_{{patch}}.log"
        )
    shell:
        """
        unset LD_LIBRARY_PATH
        
        sopa segmentation baysor {params.zarr:q} \
            --config {params.config:q} \
            --patch-index {wildcards.patch} \
            --min-area {params.min_area} \
            > {log:q} 2>&1
        """


rule sopa_resolve_baysor:
    input:
        transcripts_patch_outputs
    output:
        sentinel=touch(
            str(
                SEGMENT_DIR /
                f"patches/{{sample}}/baysorResolve_{BS_MASK_KEY}_{BS_HASH}.sentinel"
            )
        )
    params:
        zarr=segment_zarr,
        gene=config["bs_gene_column"],
        mask_key=BS_MASK_KEY,
    resources:
        **RESOURCES,
    log:
        str(
            SEGMENT_DIR /
            f"logs/{{sample}}_sopaResolve_{BS_MASK_KEY}_{BS_HASH}.log"
        )
    run:
        zarr = Path(params.zarr)
        source = zarr / "shapes" / "baysor_boundaries"
        destination = zarr / "shapes" / params.mask_key

        if source != destination and destination.exists():
            raise FileExistsError(
                f"Segmentation already exists: {destination}"
            )

        patch_dir = Path(params.zarr) / ".sopa_cache" / "transcript_patches"

        for patch in patch_dir.iterdir():
            segmentation_dir = patch / "segmentation"

            if not segmentation_dir.is_dir():
                continue

            for item in segmentation_dir.iterdir():
                destination = patch / item.name

                if destination.exists():
                    raise FileExistsError(f"Destination already exists: {destination}")

                shutil.move(str(item), str(destination))

            segmentation_dir.rmdir()
        
        command = [
            "sopa", "resolve", "baysor", str(zarr),
            "--gene-column", params.gene,
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

        #if source != destination:
        #    source.rename(destination)

        with open(log[0], "a") as handle:
            handle.write(
                f"\nRegistered segmentation: {destination}\n"
            )

