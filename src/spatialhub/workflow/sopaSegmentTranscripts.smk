
"""
Transcripts-based cell segmentation using SOPA (Baysor).

Input:
    Existing SpatialData Zarr stores.

Output:
    Reduced SpatialData Zarr stores containing segmentation results.
    => segmentation mask registered back to original Zarr store

Stages:
    1. Create reduced segmentation Zarr.
    2. Generate SOPA transcripts patches.
    3. Run Baysor segmentation per patch.
    4. Resolve segmentation boundaries (using SOPA).
    5. Register segmentation mask (shapes) to original Zarr store.
"""

from pathlib import Path
import hashlib
import json
import shlex
import subprocess
import sys

import pandas as pd
import spatialdata as spd


# ----------------------- Global configuration ----------------------- #

configfile: "sopaSegmentTranscripts.yaml"

RESOURCES = {"threads": 4, "mem_mb": 16000, "time": "04:00:00", "partition": "short"}
RESOURCES.update(config.get("resources", {}))

ZARR_DIR = Path(config["zarr_dir"])
SEGMENT_DIR = Path(config.get("segment_dir", "sopa.segment.dir"))

# Transcripts patching parameters
PATCH_WIDTH = config["tx_patch_width"]
PATCH_OVERLAP = config["tx_patch_ovlp"]


# ----------------------- Baysor configuration ----------------------- #

BS_CONFIG = config["bs_config"]
BS_MASK_KEY = config["bs_mask_key"]

BS_PARAMS = {
    "patch_width": PATCH_WIDTH,
    "patch_overlap": PATCH_OVERLAP,
    "prior": config.get("prior_segment_key"),
    "baysor_config": BS_CONFIG,
    "baysor_area": config["bs_cell_area"],
    "mask_key": config["bs_mask_key"],
}

BS_HASH = hashlib.sha256(
    json.dumps(BS_PARAMS, sort_keys=True).encode()
).hexdigest()[:12]


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

rule full:
    input:
        expand(
            str(
                SEGMENT_DIR /
                f"logs/{{sample}}_sopaRegister_baysor_{BS_HASH}.sentinel"
            ),
            sample=SAMPLE_IDS,
        )
    resources:
        **RESOURCES,
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


# ----------------------- Patch validation ----------------------- #

rule count_patch_transcripts:
    input:
        patch_indices=str(
            SEGMENT_DIR / "patches/{sample}_transcripts.indices"
        )
    output:
        counts=str(
            SEGMENT_DIR / "patches/{sample}_transcriptCounts.tsv"
        )
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
        ),
        patch_indices=str(
            SEGMENT_DIR / "patches/{sample}_transcripts.indices"
        )
    output:
        validated=str(
            SEGMENT_DIR / "patches/{sample}_transcriptCounts.validated"
        )
    resources:
        **RESOURCES,
    params:
        minimum=int(config.get("tx_min_transcripts", 4000))
    run:
        counts = pd.read_csv(input.counts, sep="\t")

        # Check whether the table exists and contains the expected data
        if counts.empty:
            raise ValueError("No transcript patch counts found.")

        if counts["n_transcripts"].isna().any():
            raise ValueError("Missing transcript counts detected.")

        expected = {
            int(p)
            for p in Path(input.patch_indices).read_text().strip().split(",")
        }
        observed = counts["patch_index"]

        if observed.isna().any() or observed.duplicated().any():
            raise ValueError("Invalid or duplicate patch indices in count table.")

        if set(observed.astype(int)) != expected:
            raise ValueError(
                "Transcript count table does not match SOPA patch indices."
            )

        # Check whether all patches have enough transcripts
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


def get_patch_indices(wildcards):
    ckpt = checkpoints.sopa_transcripts_patch.get(
        sample=wildcards.sample
    )

    with open(ckpt.output.patch_indices) as handle:
        patch_indices = [
            int(p) for p in handle.read().strip().split(",")
        ]

    return patch_indices


def baysor_patch_outputs(wildcards):
    return expand(
        str(
            SEGMENT_DIR /
            f"patches/{{sample}}/baysorSegment_{BS_MASK_KEY}_{BS_HASH}_{{patch}}.sentinel"
        ),
        sample=wildcards.sample,
        patch=get_patch_indices(wildcards),
    )


# ----------------------- Baysor segmentation ----------------------- #

# Track BS_HASH configuration to prevent later resolution with mixed patches
rule track_baysor:
    input:
        validated=str(
            SEGMENT_DIR / "patches/{sample}_transcriptCounts.validated"
        ),
    output:
        current=str(SEGMENT_DIR / "patches/{sample}_baysor_current.json")
    params:
        hash=BS_HASH
    resources:
        **RESOURCES
    run:
        Path(output.current).write_text(json.dumps({"hash": params.hash}) + "\n")


rule sopa_segment_baysor:
    input:
        current=str(SEGMENT_DIR / "patches/{sample}_baysor_current.json")
    output:
        sentinel=touch(
            str(SEGMENT_DIR / f"patches/{{sample}}/baysorSegment_{BS_MASK_KEY}_{BS_HASH}_{{patch}}.sentinel")
        )
    params:
        zarr=segment_zarr,
        config=repr(config["bs_config"]),
        min_area=config["bs_cell_area"],
        patch_dir=lambda wc: (
            Path(segment_zarr(wc))
            / ".sopa_cache"
            / "transcript_patches"
            / wc.patch
        ),
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
        : > {log:q}

        for file in {params.patch_dir:q}/segmentation* {params.patch_dir:q}/config.toml; do
            if [ -f "$file" ]; then
                echo "WARNING: Removing existing Baysor output: $file" \
                    | tee -a {log:q} >&2
                rm -- "$file"
            fi
        done

        unset LD_LIBRARY_PATH

        sopa segmentation baysor {params.zarr:q} \
            --config {params.config:q} \
            --patch-index {wildcards.patch} \
            --min-area {params.min_area} \
            >> {log:q} 2>&1
        """


def check_owner(path):
    actual = json.loads(Path(path).read_text()).get("hash")
    if actual != BS_HASH:
        raise RuntimeError(
            f"Baysor patch cache belongs to {actual}, not {BS_HASH}. "
            "Rerun segmentation with the intended configuration."
        )

rule record_completion:
    input:
        patches=baysor_patch_outputs,
        current=str(SEGMENT_DIR / "patches/{sample}_baysor_current.json")
    output:
        completed=str(
            SEGMENT_DIR /
            f"configs/{{sample}}/baysorComplete_{BS_HASH}.json"
        )
    resources:
        **RESOURCES
    run:
        check_owner(input.current)
        
        state = {
            "hash": BS_HASH,
            "status": "complete",
            "parameters": BS_PARAMS,
        }

        Path(output.completed).write_text(
            json.dumps(state, indent=2, sort_keys=True) + "\n"
        )


# ----------------------- Baysor resolution ----------------------- #

rule sopa_resolve_baysor:
    input:
        completed=str(
            SEGMENT_DIR /
            f"configs/{{sample}}/baysorComplete_{BS_HASH}.json"
        ),
        current=str(
            SEGMENT_DIR / "patches/{sample}_baysor_current.json"
        )
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
        check_owner(input.current)
        
        completion = json.loads(Path(input.completed).read_text())
        if (
            completion.get("hash") != BS_HASH
            or completion.get("status") != "complete"
            or completion.get("parameters") != BS_PARAMS
        ):
            raise ValueError("Baysor patch completion record is invalid.")
        
        zarr = Path(params.zarr)

        command = [
            "sopa", "resolve", "baysor", str(zarr),
            "--gene-column", params.gene,
            "--min-area", str(config["bs_cell_area"]),
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

rule sopa_register_mask:
    input:
        sentinel=str(
            SEGMENT_DIR /
            f"patches/{{sample}}/baysorResolve_{BS_MASK_KEY}_{BS_HASH}.sentinel"
        ),
        current=str(SEGMENT_DIR / "patches/{sample}_baysor_current.json")
    output:
        sentinel=touch(
            str(
                SEGMENT_DIR /
                f"logs/{{sample}}_sopaRegister_baysor_{BS_HASH}.sentinel"
            )
        )
    params:
        zarr=segment_zarr,
        original_zarr=input_zarr,
        source_key="baysor_boundaries",
        mask_key=BS_MASK_KEY,
        image_key=config.get("image_key", "image"),
        overwrite=config.get("overwrite_mask", False),
    resources:
        **RESOURCES
    log:
        str(
            SEGMENT_DIR /
            f"logs/{{sample}}_sopaRegister_baysor_{BS_HASH}.log"
        )
    run:
        check_owner(input.current)
        
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