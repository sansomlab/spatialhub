
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
import shutil
import subprocess
import hashlib
import json
import spatialdata as spd
import pandas as pd


# ----------------------- Configuration ----------------------- #

configfile: "sopaSegmentTranscripts.yaml"

RESOURCES = {"threads": 4, "mem_mb": 16000, "time": "04:00:00", "partition": "short"}
RESOURCES.update(config.get("resources", {}))

ZARR_DIR = Path(config["zarr_dir"])
SEGMENT_DIR = Path(config.get("segment_dir", "sopa.segment.dir"))

RUN_BAYSOR = config.get("run_baysor", False)
RUN_PROSEG = config.get("run_proseg", False)

# For beta release only!
if RUN_BAYSOR and RUN_PROSEG:
    raise ValueError(
        "Running Baysor and Proseg together is not supported "
        "in the beta release. Enable only one method."
    )

# Shared image patching parameters
PATCH_WIDTH = config["tx_patch_width"]
PATCH_OVERLAP = config["tx_patch_ovlp"]

BS_PARAMS = {
    "patch_width": PATCH_WIDTH,
    "patch_overlap": PATCH_OVERLAP,
    "baysor_config": config["bs_config"],
    "baysor_area": config["bs_cell_area"],
    "mask_key": config["bs_mask_key"],
}

BS_HASH = hashlib.sha256(
    json.dumps(BS_PARAMS, sort_keys=True).encode()
).hexdigest()[:12]

BS_MASK_KEY = config["bs_mask_key"]


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

if RUN_BAYSOR:
    REGISTRATION["baysor"] = {
        "source_key": "baysor_boundaries",
        "mask_key": BS_MASK_KEY,
        "hash": BS_HASH,
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
        **RESOURCES,
    shell:
        """
        echo "=======> All done! <======="
        """


# ----------------------- Segmentation setup ----------------------- #

#RESOURCES["mem_mb"] = config["sdata_mem"] 

rule segment_setup:
    output:
        sentinel=touch(
            str(SEGMENT_DIR / "logs/{sample}_sopa_segmentSetup.sentinel")
        )
    params:
        original_zarr=input_zarr,
        zarr=segment_zarr,
        mode="transcripts",
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


# ----------------------- Transcripts patching ----------------------- #

checkpoint sopa_transcripts_patch:
    input:
        str(SEGMENT_DIR / "logs/{sample}_sopa_segmentSetup.sentinel")
    output:
        patch_indices=str(
            SEGMENT_DIR / "patches/{sample}_transcripts.indices"
        )
    params:
        zarr=segment_zarr,
        width=PATCH_WIDTH,
        overlap=PATCH_OVERLAP,
        prior=config.get("prior_segment_key")
    resources:
        **RESOURCES,
    log:
        str(SEGMENT_DIR / "logs/{sample}_patchTranscripts.log")
    run:
        command = [
            "sopa", "patchify", "transcripts", params.zarr,
            "--patch-width-microns", str(params.width),
            "--patch-overlap-microns", str(params.overlap),
        ]
        if params.prior:
            command.extend([
                "--prior-shapes-key", str(params.prior)
            ])

        with open(log[0], "w") as handle:
            subprocess.run(
                command,
                stdout=handle,
                stderr=subprocess.STDOUT,
                check=True,
            )

        # Retrieve actual transcript patch indices
        sdata = spd.read_zarr(params.zarr)
        patch_indices = sdata["transcripts_patches"].index.tolist()

        if not patch_indices:
            raise ValueError("No transcript patches were generated.")

        with open(output.patch_indices, "w") as handle:
            handle.write(",".join(map(str, patch_indices)))


rule count_patch_transcripts:
    input:
        patch_indices=str(
            SEGMENT_DIR / "patches/{sample}_transcripts.indices"
        )
    output:
        counts=str(
            SEGMENT_DIR / "patches/{sample}_transcriptCounts.tsv"
        )
    log:
        str(SEGMENT_DIR / "logs/{sample}_countPatchTx.log")
    params:
        zarr=segment_zarr,
    threads: 4
    resources:
        **RESOURCES,
    log:
        str(SEGMENT_DIR / "logs/{sample}_transcriptCounts.log")
    shell:
        """
        python -m spatialhub.scripts.countPatchTranscripts \
            --input-zarr {params.zarr:q} \
            --output {output.counts:q} \
            --workers {threads} \
            > {log:q} 2>&1
        """


rule validate_patch_transcripts:
    input:
        counts=str(
            SEGMENT_DIR / "patches/{sample}_transcriptCounts.tsv"
        )
    output:
        validated=str(
            SEGMENT_DIR / "patches/{sample}_transcriptCounts.validated"
        )
    log:
        str(SEGMENT_DIR / "logs/{sample}_validatePatchTx.log")
    resources:
        **RESOURCES,
    params:
        minimum=int(config.get("tx_min_transcripts", 4000))
    run:
        counts = pd.read_csv(input.counts, sep="\t")

        if counts.empty:
            raise ValueError("No transcript patch counts found.")

        if counts["n_transcripts"].isna().any():
            raise ValueError("Missing transcript counts detected.")

        insufficient = counts[
            counts["n_transcripts"] < params.minimum
        ]

        if not insufficient.empty:
            details = ", ".join(
                f"{row.patch_index} ({row.n_transcripts})"
                for row in insufficient.itertuples()
            )

            raise ValueError(
                f"{len(insufficient)} patches contain fewer than "
                f"{params.minimum} transcripts: {details}. "
                "Increase tx_patch_width and regenerate patches."
            )

        Path(output.validated).write_text(
            f"All patches contain at least {params.minimum} transcripts.\n"
        )


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
        patch_indices=str(
            SEGMENT_DIR / "patches/{sample}_transcripts.indices"
        ),
        validated=str(
            SEGMENT_DIR / "patches/{sample}_transcriptCounts.validated"
        ),
    output:
        sentinel=touch(
            str(SEGMENT_DIR / f"patches/{{sample}}/baysorSegment_{BS_MASK_KEY}_{BS_HASH}_{{patch}}.sentinel")
        )
    params:
        zarr=segment_zarr,
        config=repr(config["bs_config"]),
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
    resources:
        **RESOURCES
    log:
        str(
            SEGMENT_DIR /
            f"logs/{{sample}}_sopaResolve_{BS_MASK_KEY}_{BS_HASH}.log"
        )
    run:
        zarr = Path(params.zarr)
        patch_dir = zarr / ".sopa_cache" / "transcript_patches"

        # Normalize Baysor 0.9.0 output locations
        # (disabled for Baysor v0.7.1)
        #for patch in patch_dir.iterdir():
        #    segmentation_dir = patch / "segmentation"
        #
        #    if not segmentation_dir.is_dir():
        #        continue
        #
        #    for item in segmentation_dir.iterdir():
        #        destination = patch / item.name
        #
        #        if destination.exists():
        #            raise FileExistsError(
        #                f"Destination already exists: {destination}"
        #            )
        #
        #        shutil.move(str(item), str(destination))
        #
        #    segmentation_dir.rmdir()

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

        # Verify that SOPA produced the resolved shape.
        sdata = spd.read_zarr(zarr)

        if "baysor_boundaries" not in sdata.shapes:
            raise RuntimeError(
                "Baysor resolution did not produce "
                "'baysor_boundaries'."
            )


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
        method="baysor|proseg",
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