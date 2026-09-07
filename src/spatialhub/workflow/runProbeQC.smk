"""
Probe QC workflow
-----------------

This workflow enables calculating quality control metrics to assess the quality of the run from the transcriptomics side.
(This workflow does *not* readily support QC of imaging or segmentation.)

Some QC metrics may be technology specific - e.g. FOV level QC metrics for CosMx. 


Inputs
------
- cell x gene transcripts count table derived from a given segmentation mask


Outputs
-------
- A CSV table of QC metrics
"""

import os

from importlib.resources import files
from pathlib import Path

import anndata as ad


configfile: "runProbeQC.yaml"


# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------

# [NOTE] this key is set by the CLI, not the config file!
if "lock" in config:
    LOCK = config["lock"].lower()
else:
    LOCK = "false"

H5AD_DIR = Path(config["h5ad_dir"])
PROBEQC_DIR = Path(config.get("probeqc_dir", "probeqc.dir"))
RUN_NAME_OVERRIDE = config.get("run_name")

RESOURCES = {
    "threads": 4,
    "mem_mb": 16000,
    "time": "04:00:00",
    "partition": "short",
}
RESOURCES.update(config.get("resources", {}))

# R script distributed with the spatialhub package
PROBEQC_SCRIPT = files("spatialhub").joinpath(
    "scripts",
    "calculateQC_fromCountsMatrix.R",
)

if config.get("runBrukerFOVqc", False):
    if not config.get("bruker_code") or not config.get("bruker_data"):
        raise ValueError(
            "runBrukerFOVqc is enabled but bruker_code and/or bruker_data "
            "are not specified."
        )


# ----------------------------------------------------------------------
# Discover input files
# ----------------------------------------------------------------------

if not H5AD_DIR.exists():
    raise FileNotFoundError(
        f"Input h5ad_dir does not exist: {H5AD_DIR}"
    )
if not H5AD_DIR.is_dir():
    raise NotADirectoryError(
        f"Expected h5ad_dir to be a directory: {H5AD_DIR}"
    )

H5AD_FILES = sorted(H5AD_DIR.glob("*.h5ad"))
RDS_FILES = sorted(H5AD_DIR.glob("*.rds"))
if H5AD_FILES and RDS_FILES:
    raise ValueError(
        "h5ad_dir contains both .h5ad and .rds files. "
        "Please use a directory containing only one input format."
    )

if H5AD_FILES:
    INPUT_FILES = H5AD_FILES
    INPUT_FORMAT = "h5ad"
elif RDS_FILES:
    INPUT_FILES = RDS_FILES
    INPUT_FORMAT = "rds"
else:
    raise FileNotFoundError(
        f"No .h5ad or .rds files found in {H5AD_DIR}"
    )

# Each filename stem defines the sample name used by Snakemake.
SAMPLES = [path.stem for path in INPUT_FILES]

if len(SAMPLES) != len(set(SAMPLES)):
    raise ValueError(
        "Duplicate sample names found after removing file extensions."
    )

INPUT_BY_SAMPLE = {
    path.stem: str(path)
    for path in INPUT_FILES
}


# ----------------------------------------------------------------------
# Validate inputs and resolve output directory name
# ----------------------------------------------------------------------

def format_points_from(value):
    """Return points_from as a stable tuple of strings."""

    if isinstance(value, str):
        return (value,)

    return tuple(str(x) for x in value)


def read_counts_mtx_config(path):
    """
    Read the aggregation configuration (i.e. source Points elements, segementation mask, etc.)
    from an H5AD AnnData table generated with `spatialhub`.
    """

    adata = ad.read_h5ad(path, backed="r")

    try:
        if "counts_mtx_source" not in adata.uns:
            return None

        source = adata.uns["counts_mtx_source"]

        required = [
            "points_from",
            "shapes_by",
            "agg_func",
            "coords_system",
        ]

        missing = [
            key
            for key in required
            if key not in source
        ]

        if missing:
            raise ValueError(
                f"{path} has incomplete "
                "adata.uns['counts_mtx_source'] metadata. "
                f"Missing: {missing}"
            )

        return {
            "points_from": format_points_from(
                source["points_from"]
            ),
            "shapes_by": str(source["shapes_by"]),
            "agg_func": str(source["agg_func"]),
            "coords_system": str(source["coords_system"]),
        }

    finally:
        adata.file.close()


def validate_h5ad_provenance(h5ad_files):
    """
    Validate counts_mtx_source metadata across all H5AD inputs.

    Returns the shared count-matrix configuration if provenance is
    available, otherwise None.

    If some files contain provenance and others do not, fail.
    """

    signatures = [
        read_counts_mtx_config(path)
        for path in h5ad_files
    ]

    have_provenance = [
        signature is not None
        for signature in signatures
    ]

    # Mixed H5AD sources
    if any(have_provenance) and not all(have_provenance):
        missing = [
            str(path)
            for path, signature in zip(h5ad_files, signatures)
            if signature is None
        ]

        raise ValueError(
            "Some H5ADs contain adata.uns['counts_mtx_source'] "
            "metadata and others do not.\n"
            "Missing provenance:\n"
            + "\n".join(missing)
        )

    # External H5AD collection (not generated using `spatialhub`)
    if not any(have_provenance):
        return None

    reference = signatures[0]

    for path, signature in zip(h5ad_files[1:], signatures[1:]):
        if signature != reference:
            raise ValueError(
                "Input H5ADs do not share the same "
                "count-matrix configuration.\n"
                f"Reference: {reference}\n"
                f"{path}: {signature}"
            )

    return reference


def set_run_name(mtx_provenance):
    """Construct the default ProbeQC run name from H5AD provenance."""

    points_label = "+".join(mtx_provenance["points_from"])

    return ".".join(
        [
            points_label,
            mtx_provenance["shapes_by"],
            mtx_provenance["agg_func"],
            mtx_provenance["coords_system"],
        ]
    )


# Apply functions to validate input and set outdir run name
if INPUT_FORMAT == "h5ad":

    PROVENANCE = validate_h5ad_provenance(H5AD_FILES)

    if RUN_NAME_OVERRIDE:
        # run_name affects only the output directory.
        RUN_NAME = str(RUN_NAME_OVERRIDE)

    elif PROVENANCE is not None:
        RUN_NAME = set_run_name(PROVENANCE)

    else:
        raise ValueError(
            "Input H5ADs do not contain "
            "adata.uns['counts_mtx_source'] metadata. "
            "Specify 'run_name' in probeQC.yaml."
        )

else:
    # RDS files do not provide the AnnData provenance mechanism,
    # so the run directory must be named explicitly.
    if not RUN_NAME_OVERRIDE:
        raise ValueError(
            "run_name must be specified in probeQC.yaml "
            "when using RDS input."
        )

    RUN_NAME = str(RUN_NAME_OVERRIDE)


RUN_DIR = PROBEQC_DIR / RUN_NAME

print(f"ProbeQC input directory: {H5AD_DIR}")
print(f"Input format: {INPUT_FORMAT}")
print(f"Samples found: {len(SAMPLES)}")
print(f"ProbeQC run name: {RUN_NAME}")
print(f"ProbeQC output directory: {RUN_DIR}")


# ----------------------------------------------------------------------
# Define output
# ----------------------------------------------------------------------

RUN_BRUKER_FOV_QC = config.get("runBrukerFOVqc", False)

PROBEQC_OUTPUTS = [
    os.path.join(
        str(RUN_DIR),
        "{sample}_cellQCmetrics.csv",
    ),
]

if RUN_BRUKER_FOV_QC:
    PROBEQC_OUTPUTS.extend(
        [
            os.path.join(
                str(RUN_DIR),
                "{sample}_fovQCmetrics.csv",
            ),
            os.path.join(
                str(RUN_DIR),
                "{sample}_instrGeneBias.csv",
            ),
        ]
    )



# ----------------------------------------------------------------------
# Workflow
# ----------------------------------------------------------------------

rule full:
    input:
        expand(
            PROBEQC_OUTPUTS,
            sample=SAMPLES,
        )


rule probeQC:
    input:
        lambda wc: INPUT_BY_SAMPLE[wc.sample],
    output:
        PROBEQC_OUTPUTS,
    log:
        os.path.join(
            str(RUN_DIR),
            "logs",
            "{sample}.probeQC.log",
        ),
    resources:
        **RESOURCES,
    params:
        outdir=str(RUN_DIR),
        background=config.get("background", "Negative"),
        poscontrol=config.get("poscontrol", "System"),
        custom_panel=lambda wc: (
            f'--customPanel "{config["custom_panel_file"]}"'
            if config.get("custom_panel_file")
            else ""
        ),
        bruker_qc=lambda wc: (
            f'--runBrukerFOVqc TRUE '
            f'--brukerCodeFile "{config["bruker_code"]}" '
            f'--brukerDataFile "{config["bruker_data"]}"'
            if config.get("runBrukerFOVqc", False)
            else ""
        ),
    shell:
        """
        mkdir -p "{RUN_DIR}" "$(dirname "{log}")"

        Rscript "{PROBEQC_SCRIPT}" \
            --path2sce "{input}" \
            --outdir "{RUN_DIR}" \
            --background "{params.background}" \
            --poscontrol "{params.poscontrol}" \
            {params.custom_panel} \
            {params.bruker_qc} \
            >"{log}" 2>&1
        """
