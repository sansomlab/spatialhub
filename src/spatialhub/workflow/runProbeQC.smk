"""
Probe QC workflow
-----------------

This workflow enables calculating metrics to assess the quality of the run from the transcriptomics side.
(This workflow does *not* readily support QC of imaging or segmentation.)

Some QC metrics may be technology-specific - e.g. FOV level QC metrics for CosMx. 


Inputs
------
- instance x gene transcripts count table derived from a given segmentation mask


Outputs
-------
- CSV tables of QC metrics and flags per sample
- Enriched instance-level metadata table for the project, including QC metrics and sample-level metadata relevant for the study
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
    "probeqc_instance_report.Rmd",
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
if not H5AD_DIR.is_dir():
    raise NotADirectoryError(
        f"Invalid H5AD directory: {H5AD_DIR}"
    )

H5AD_FILES = sorted(H5AD_DIR.glob("*.h5ad"))
if not H5AD_FILES:
    raise ValueError(f"No H5AD files found in {H5AD_DIR}")

# Each filename stem defines the sample name used by Snakemake.
SAMPLES = [path.stem for path in H5AD_FILES]


# ----------------------------------------------------------------------
# Validate inputs and resolve output directory name
# ----------------------------------------------------------------------

def read_provenance(path):
    """Read spatialhub provenance without loading the count matrix."""
    
    adata = ad.read_h5ad(path, backed="r")
    source = adata.uns.get("counts_mtx_source")
    adata.file.close()
    
    if source is not None:
        source = {
            k: v.tolist() if hasattr(v, "tolist") else v
            for k, v in source.items()
        }
    
    return source

sources = [read_provenance(path) for path in H5AD_FILES]

if any(s is not None for s in sources) and not all(s is not None for s in sources):
    raise ValueError("Cannot mix spatialhub and external H5AD inputs")

SPATIALHUB = all(s is not None for s in sources)

if SPATIALHUB:
    # Sample-specific provenance is expected to differ. Other fields must match.
    def signature(source):
        return {
            k: v for k, v in source.items()
            if k not in ("sample_id", "input_zarr")
        }

    PROVENANCE = signature(sources[0])
    for path, source in zip(H5AD_FILES[1:], sources[1:]):
        if signature(source) != PROVENANCE:
            raise ValueError(f"Incompatible count-matrix provenance: {path}")

    region_type = PROVENANCE.get("region_type")
    if region_type == "segment":
        region_name = PROVENANCE["shapes_by"]
    elif region_type == "grid":
        region_name = f"{PROVENANCE['bin_size']}{PROVENANCE['bin_unit']}"
    else:
        raise ValueError(f"Unknown region_type: {region_type}")

    DEFAULT_RUN_NAME = f"{region_type}_{region_name}"
else:
    DEFAULT_RUN_NAME = None


# Set run name and output directories accordingly
RUN_NAME = config.get("run_name", DEFAULT_RUN_NAME)

if not RUN_NAME:
    raise ValueError(
        "H5AD inputs generated outside of spatialhub require an explicit run_name"
    )

RUN_DIR = Path(config["outdir"]) / RUN_NAME
RES_DIR = RUN_DIR / "results"
QC_DIR = RUN_DIR / "filtering" / "current"
REPORT_DIR = QC_DIR / "reports"

print(f"ProbeQC input directory: {H5AD_DIR}")
print(f"Input source: {'spatialhub' if SPATIALHUB else 'external'}")
print(f"Samples found: {len(SAMPLES)}")
print(f"ProbeQC output directory: {RUN_DIR}")


if RUN_BRUKER_QC and SPATIALHUB and region_type == "grid":
    raise ValueError(
        "Bruker FOV QC is not supported for grid-based count matrices. "
        "Set run_bruker_qc: false."
    )


# ----------------------------------------------------------------------
# Define output
# ----------------------------------------------------------------------

MTXQC_OUTPUTS = [
    os.path.join(
        str(RES_DIR),
        "{sample}_instanceQCmetrics.csv",
    ),
    os.path.join(
        str(RES_DIR),
        "{sample}_instanceMetadata.csv",
    ),
    os.path.join(
        str(RES_DIR),
        "{sample}_probeQCmetrics.csv",
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
    os.path.join(str(QC_DIR), f"{RUN_NAME}_instanceQCmetrics.csv"),
    os.path.join(str(QC_DIR), f"{RUN_NAME}_sampleQCmetrics.csv")
]

if RUN_FOV_QC:
    FLAG_OUTPUTS.append(
        os.path.join(str(QC_DIR), f"{RUN_NAME}_fovQCmetrics.csv")
    )

REPORT_OUTPUTS = [
    os.path.join(str(REPORT_DIR), f"probeqc_{RUN_NAME}_instance.html"),
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
        results_dir=str(RES_DIR)
    shell:
        """
        if [ {params.lock} = true ]; then
            chmod -R a-w "{params.results_dir}"
        fi
        echo "=======> All done! <======="
        """


rule probeQC_mtx:
    input:
        str(H5AD_DIR / "{sample}.h5ad"),
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
            f"{RUN_NAME}_instanceQCmetrics.csv",
        ),
        thresholds=str(QC_THRESHOLDS),
        report_source=str(PROBEQC_REPORT_RMD),
        header_source=str(PROBEQC_HEADER_RMD),
        graphical_params_source=str(PROBEQC_GRAPHICAL_PARAMS_RMD),
    output:
        html=os.path.join(
            str(REPORT_DIR),
            f"probeqc_{RUN_NAME}_instance.html",
        ),
    params:
        report_dir=str(REPORT_DIR),
        run_dir=str(QC_DIR),
        group=config["group"],

        # QC thresholds
        instance_nCount_cutoff=qc_thresholds["instance_nCount_cutoff"],
        instance_nFeature_cutoff=qc_thresholds["instance_nFeature_cutoff"],
        instance_percentNeg_cutoff=qc_thresholds["instance_percentNeg_cutoff"],
        instance_complexity_cutoff=qc_thresholds["instance_complexity_cutoff"],
        sampleMinInstances=qc_thresholds["sampleMinInstances"],
        sampleMinPercent=qc_thresholds["sampleMinPercent"],
        exclInstancesInFailedFOVs=(
            "TRUE"
            if qc_thresholds["exclInstancesInFailedFOVs"]
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

        cp --update "{input.report_source}" "{params.run_dir}/probeqc_instance_report.Rmd"
        cp --update "{input.header_source}" "{params.run_dir}/probeqc_header.Rmd"
        cp --update "{input.graphical_params_source}" "{params.run_dir}/_graphical_params.Rmd"

        Rscript -e '
        rmarkdown::render(
            input = "{params.run_dir}/probeqc_instance_report.Rmd",
            output_file = "probeqc_{RUN_NAME}_instance.html",
            output_dir = "{params.report_dir}",
            params = list(
                group = "{params.group}",
                instance_nCount_cutoff = {params.instance_nCount_cutoff},
                instance_nFeature_cutoff = {params.instance_nFeature_cutoff},
                instance_percentNeg_cutoff = {params.instance_percentNeg_cutoff},
                instance_complexity_cutoff = {params.instance_complexity_cutoff},
                sampleMinInstances = {params.sampleMinInstances},
                sampleMinPercent = {params.sampleMinPercent},
                exclInstancesInFailedFOVs = {params.exclInstancesInFailedFOVs},
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
            f"{RUN_NAME}_instanceQCmetrics.csv",
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

        cp --no-clobber "{input.report_source}" "{params.run_dir}/probeqc_fov_report.Rmd"
        cp --no-clobber "{input.bruker_source}" "{params.run_dir}/probeqc_bruker.Rmd"
        cp --no-clobber "{input.header_source}" "{params.run_dir}/probeqc_header.Rmd"
        cp --no-clobber "{input.graphical_params_source}" "{params.run_dir}/_graphical_params.Rmd"

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