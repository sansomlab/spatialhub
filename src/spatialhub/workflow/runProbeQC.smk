"""
Probe QC workflow
-----------------

This workflow enables calculating metrics to assess the quality of the run from the transcriptomics side.
(This workflow does *not* readily support QC of imaging or segmentation.)

Some QC metrics may be technology-specific - e.g. FOV level QC metrics for CosMx. 


Inputs
------
- cell x gene transcripts count table derived from a given segmentation mask


Outputs
-------
- CSV tables of QC metrics and flags per sample
- Enriched cell-level metadata table for the project, including QC metrics and sample-level metadata relevant for the study
- HTML QC reports and underlying Rmd files for the user to further amend manually
"""

import os
from importlib.resources import files
from pathlib import Path
import yaml
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

RESOURCES = {
    "threads": 4,
    "mem_mb": 16000,
    "time": "04:00:00",
    "partition": "short",
}
RESOURCES.update(config.get("resources", {}))


# Pointers to R scripts distributed with the spatialhub package
PROBEQC_MTX_SCRIPT = files("spatialhub").joinpath(
    "scripts",
    "calculateQC_fromCountsMatrix.R",
)
PROBEQC_DERIVED_SCRIPT = files("spatialhub").joinpath(
    "scripts",
    "calculateQC_derivedVariables.R",
)
PROBEQC_FLAG_SCRIPT = files("spatialhub").joinpath(
    "scripts",
    "flagProbeQCfailures.R",
)

# Pointers to Rmd files distributed within the spatialhub package
PROBEQC_REPORT_RMD = files("spatialhub").joinpath(
    "reports",
    "probeqc_cell_report.Rmd",
)
PROBEQC_REPORT_FOV_RMD = files("spatialhub").joinpath(
    "reports",
    "probeqc_fov_report.Rmd",
)
PROBEQC_BRUKER_RMD = files("spatialhub").joinpath(
    "reports",
    "probeqc_bruker.Rmd",
)
PROBEQC_HEADER_RMD = files("spatialhub").joinpath(
    "reports",
    "probeqc_header.Rmd",
)
PROBEQC_GRAPHICAL_PARAMS_RMD = files("spatialhub").joinpath(
    "reports",
    "_graphical_params.Rmd",
)


# Pre-flight checks for potential missing information
RUN_BRUKER_QC = config.get("run_bruker_qc", False)
RUN_FOV_QC = config.get("run_fov_qc", False)

if RUN_BRUKER_QC and not RUN_FOV_QC:
    raise ValueError(
        "run_bruker_qc requires run_fov_qc to be enabled."
    )

if RUN_BRUKER_QC:
    if not config.get("bruker_code") or not config.get("bruker_data"):
        raise ValueError(
            "run_bruker_qc is enabled but bruker_code and/or bruker_data "
            "are not specified."
        )


# Fetching QC thresholds from mini YAML
QC_THRESHOLDS = Path(config["qc_thresh"])

with QC_THRESHOLDS.open() as f:
    qc_thresholds = yaml.safe_load(f)




# ----------------------------------------------------------------------
# Discover input files
# ----------------------------------------------------------------------

H5AD_DIR = Path(config["h5ad_dir"])

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

# Initiate global output directory names - with default if not specified
PROBEQC_DIR = Path(config.get("probeqc_dir", "probeqc.dir"))
RUN_NAME_OVERRIDE = config.get("run_name")
QC_PROFILE = config.get("qc_outdir", "current")


def format_points_from(value):
    """Return points_from as a stable tuple of strings."""

    if isinstance(value, str):
        return (value,)

    return tuple(str(x) for x in value)


def read_counts_mtx_config(path):
    """
    Read the aggregation configuration (i.e. source Points elements, segmentation mask, etc.)
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
        # `run_name` affects only the output *directory* name
        RUN_NAME = str(RUN_NAME_OVERRIDE)

    elif PROVENANCE is not None:
        RUN_NAME = set_run_name(PROVENANCE)

    else:
        raise ValueError(
            "Input H5ADs do not contain "
            "adata.uns['counts_mtx_source'] metadata. "
            "Specify 'run_name' in runProbeQC.yaml."
        )

else:
    # RDS files do not provide the AnnData provenance mechanism,
    # so the run directory must be named explicitly.
    if not RUN_NAME_OVERRIDE:
        raise ValueError(
            "run_name must be specified in runProbeQC.yaml "
            "when using RDS input."
        )

    RUN_NAME = str(RUN_NAME_OVERRIDE)


# Define sub-directory names based on input
RUN_DIR = PROBEQC_DIR / RUN_NAME
RES_DIR = RUN_DIR / "results" 
QC_DIR = RUN_DIR / "filtering" / QC_PROFILE
REPORT_DIR = QC_DIR / "reports"

print(f"ProbeQC input directory: {H5AD_DIR}")
print(f"Input format: {INPUT_FORMAT}")
print(f"Samples found: {len(SAMPLES)}")
print(f"ProbeQC run name: {RUN_NAME}")
print(f"ProbeQC metrics calculation directory: {RES_DIR}")
print(f"ProbeQC thresholding results directory: {QC_DIR}")


# ----------------------------------------------------------------------
# Define output
# ----------------------------------------------------------------------

MTXQC_OUTPUTS = [
    os.path.join(
        str(RES_DIR),
        "{sample}_cellQCmetrics.csv",
    ),
    os.path.join(
        str(RES_DIR),
        "{sample}_cellMetadata.csv",
    ),
    os.path.join(
        str(RES_DIR),
        "{sample}_probeClassifier.csv",
    ),
]

if RUN_BRUKER_QC:
    MTXQC_OUTPUTS.extend(
        [
            os.path.join(
                str(RES_DIR),
                "{sample}_brukerQCresults.csv",
            ),
        ]
    )


DERIVEDQC_OUTPUTS = [
    os.path.join(
        str(RES_DIR),
        "{sample}_sampleQCmetrics.csv",
    ),
]

if RUN_FOV_QC:
    DERIVEDQC_OUTPUTS.extend(
        [
           os.path.join(
                str(RES_DIR),
                "{sample}_fovQCmetrics.csv",
            ) 
        ]
    )


FLAG_OUTPUTS = [
    os.path.join(str(QC_DIR), f"{RUN_NAME}_cellQCmetrics.csv"),
    os.path.join(str(QC_DIR), f"{RUN_NAME}_sampleQCmetrics.csv")
]

if RUN_FOV_QC:
    FLAG_OUTPUTS.append(
        os.path.join(str(QC_DIR), f"{RUN_NAME}_fovQCmetrics.csv")
    )

REPORT_OUTPUTS = [
    os.path.join(str(REPORT_DIR), f"probeqc_{RUN_NAME}_cell.html"),
]

if RUN_FOV_QC:
    REPORT_OUTPUTS.append(
        os.path.join(str(REPORT_DIR), f"probeqc_{RUN_NAME}_fov.html")
    )




# ----------------------------------------------------------------------
# Workflow
# ----------------------------------------------------------------------

rule full:
    input:
        expand(MTXQC_OUTPUTS + DERIVEDQC_OUTPUTS, sample=SAMPLES,),
        FLAG_OUTPUTS, REPORT_OUTPUTS,
    resources:
        **RESOURCES,
    params:
        lock=LOCK,
        probeqcdir=RUN_DIR,
    shell:
        """
        if [ {params.lock} = true ]; then
            chmod -R a-w {params.probeqcdir}
        fi
        echo "=======> All done! <======="
        """


rule probeQC_mtx:
    input:
        lambda wc: INPUT_BY_SAMPLE[wc.sample],
    output:
        MTXQC_OUTPUTS,
    log:
        os.path.join(
            str(RUN_DIR),
            "logs",
            "{sample}.probeQC_mtx.log",
        ),
    resources:
        **RESOURCES,
    params:
        outdir=str(RES_DIR),
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
            if RUN_BRUKER_QC
            else ""
        ),
    shell:
        """
        mkdir -p "{RES_DIR}" "$(dirname "{log}")"

        Rscript "{PROBEQC_MTX_SCRIPT}" \
            --path2sce "{input}" \
            --outdir "{RES_DIR}" \
            --background "{params.background}" \
            --poscontrol "{params.poscontrol}" \
            {params.custom_panel} \
            {params.bruker_qc} \
            >"{log}" 2>&1
        """


rule probeQC_derived:
    input:
        lambda wc: expand(
            MTXQC_OUTPUTS,
            sample=[wc.sample],
        ),
    output:
        DERIVEDQC_OUTPUTS,
    log:
       os.path.join(
            str(RUN_DIR),
            "logs",
            "{sample}.probeQC_derived.log",
        ),
    resources:
        **RESOURCES,
    params:
        sample_id=lambda wc: wc.sample,
        run_dir=str(RES_DIR),
        run_fov_qc=lambda wc: (
            "TRUE"
            if RUN_FOV_QC
            else "FALSE"
        )
    shell:
        """
        mkdir -p "$(dirname "{log}")"

        Rscript "{PROBEQC_DERIVED_SCRIPT}" \
            --sample_id "{params.sample_id}" \
            --run_dir "{params.run_dir}" \
            --runFOVqc "{params.run_fov_qc}" \
            >"{log}" 2>&1
        """


rule flag_probeQC:
    input:
        thresholds=str(QC_THRESHOLDS),
        metrics=expand(
            MTXQC_OUTPUTS + DERIVEDQC_OUTPUTS,
            sample=SAMPLES,
        ),
    output:
        FLAG_OUTPUTS,
    log:
        os.path.join(str(QC_DIR), "probeQC_flag.log"),
    resources:
        **RESOURCES,
    params:
        run_dir=str(RES_DIR),
        metadata_file=str(config['metadata_file']),
        run_name=str(RUN_NAME),
        outdir=str(QC_DIR),
        runFOVqc=RUN_FOV_QC,
        runBrukerFOVqc=RUN_BRUKER_QC,
    shell:
        """
        mkdir -p "{params.outdir}"

        Rscript "{PROBEQC_FLAG_SCRIPT}" \
            --run_dir "{params.run_dir}" \
            --metadata_file "{params.metadata_file}" \
            --run_name "{params.run_name}" \
            --outdir "{params.outdir}" \
            --qc_thresholds "{input.thresholds}" \
            --runFOVqc "{params.runFOVqc}" \
            --runBrukerFOVqc "{params.runBrukerFOVqc}" \
            >"{log}" 2>&1
        """


rule probeQC_report:
    input:
        metrics=os.path.join(
            str(QC_DIR),
            f"{RUN_NAME}_cellQCmetrics.csv",
        ),
        thresholds=str(QC_THRESHOLDS),
        report_source=str(PROBEQC_REPORT_RMD),
        header_source=str(PROBEQC_HEADER_RMD),
        graphical_params_source=str(PROBEQC_GRAPHICAL_PARAMS_RMD),
    output:
        html=os.path.join(
            str(REPORT_DIR),
            f"probeqc_{RUN_NAME}_cell.html",
        ),
    params:
        report_dir=str(REPORT_DIR),
        run_dir=str(QC_DIR),
        group=config["group"],

        # QC thresholds
        cell_nCount_cutoff=qc_thresholds["cell_nCount_cutoff"],
        cell_nFeature_cutoff=qc_thresholds["cell_nFeature_cutoff"],
        cell_percentNeg_cutoff=qc_thresholds["cell_percentNeg_cutoff"],
        cell_complexity_cutoff=qc_thresholds["cell_complexity_cutoff"],
        sampleMinCells=qc_thresholds["sampleMinCells"],
        sampleMinPercent=qc_thresholds["sampleMinPercent"],
        exclCellsInFailedFOVs=(
            "TRUE"
            if qc_thresholds["exclCellsInFailedFOVs"]
            else "FALSE"
        ),

        # Plotting parameters
        composite_ncol=config.get("composite_ncol", 2),
        composite_spacing=config.get("composite_spacing", 500),
        plot_ncol_spatial=config.get("plot_ncol_spatial", 4),
        plot_ncol=config.get("plot_ncol", 2),
        font_size=config.get("font_size", 10),
        legend_size=config.get("legend_size", 3),

    log:
        os.path.join(
            str(QC_DIR),
            "probeQC_report.log",
        ),
    resources:
        **RESOURCES,
    shell:
        """
        mkdir -p "{params.report_dir}"

        cp --update "{input.report_source}" "{params.run_dir}/probeqc_cell_report.Rmd"
        cp --update "{input.header_source}" "{params.run_dir}/probeqc_header.Rmd"
        cp --update "{input.graphical_params_source}" "{params.run_dir}/_graphical_params.Rmd"

        Rscript -e '
        rmarkdown::render(
            input = "{params.run_dir}/probeqc_cell_report.Rmd",
            output_file = "probeqc_{RUN_NAME}_cell.html",
            output_dir = "{params.report_dir}",
            params = list(
                group = "{params.group}",
                cell_nCount_cutoff = {params.cell_nCount_cutoff},
                cell_nFeature_cutoff = {params.cell_nFeature_cutoff},
                cell_percentNeg_cutoff = {params.cell_percentNeg_cutoff},
                cell_complexity_cutoff = {params.cell_complexity_cutoff},
                sampleMinCells = {params.sampleMinCells},
                sampleMinPercent = {params.sampleMinPercent},
                exclCellsInFailedFOVs = {params.exclCellsInFailedFOVs},
                composite_ncol = {params.composite_ncol},
                composite_spacing = {params.composite_spacing},
                plot_ncol_spatial = {params.plot_ncol_spatial},
                plot_ncol = {params.plot_ncol},
                font_size = {params.font_size},
                legend_size = {params.legend_size}
            )
        )
        ' >"{log}" 2>&1
        """

rule probeQC_report_FOV:
    input:
        metrics=os.path.join(
            str(QC_DIR),
            f"{RUN_NAME}_cellQCmetrics.csv",
        ),
        thresholds=str(QC_THRESHOLDS),
        report_source=str(PROBEQC_REPORT_FOV_RMD),
        header_source=str(PROBEQC_HEADER_RMD),
        bruker_source=str(PROBEQC_BRUKER_RMD),
        graphical_params_source=str(PROBEQC_GRAPHICAL_PARAMS_RMD),
    output:
        html=os.path.join(
            str(REPORT_DIR),
            f"probeqc_{RUN_NAME}_fov.html",
        ),
    params:
        report_dir=str(REPORT_DIR),
        run_dir=str(QC_DIR),
        group=config["group"],

        # QC parameters
        run_bruker_qc=(
            "TRUE"
            if RUN_BRUKER_QC
            else "FALSE"
        ),

        # Plotting parameters
        composite_ncol=config.get("composite_ncol", 2),
        composite_spacing=config.get("composite_spacing", 500),
        plot_ncol_spatial=config.get("plot_ncol_spatial", 4),
        plot_ncol=config.get("plot_ncol", 2),
        font_size=config.get("font_size", 10),
        legend_size=config.get("legend_size", 3),

    log:
        os.path.join(
            str(QC_DIR),
            "probeQC_report_fov.log",
        ),
    resources:
        **RESOURCES,
    shell:
        """
        mkdir -p "{params.report_dir}"

        cp --update "{input.report_source}" "{params.run_dir}/probeqc_fov_report.Rmd"
        cp --update "{input.bruker_source}" "{params.run_dir}/probeqc_bruker.Rmd"
        cp --update "{input.header_source}" "{params.run_dir}/probeqc_header.Rmd"
        cp --update "{input.graphical_params_source}" "{params.run_dir}/_graphical_params.Rmd"

        Rscript -e '
        rmarkdown::render(
            input = "{params.run_dir}/probeqc_fov_report.Rmd",
            output_file = "probeqc_{RUN_NAME}_fov.html",
            output_dir = "{params.report_dir}",
            params = list(
                group = "{params.group}",
                run_bruker_qc = {params.run_bruker_qc},
                composite_ncol = {params.composite_ncol},
                composite_spacing = {params.composite_spacing},
                plot_ncol_spatial = {params.plot_ncol_spatial},
                plot_ncol = {params.plot_ncol},
                font_size = {params.font_size},
                legend_size = {params.legend_size}
            )
        )
        ' >"{log}" 2>&1
        """