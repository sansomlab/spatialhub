import os
import pandas as pd
import re


configfile: "cosmx_makeImage.yaml"


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

OUTDIR = config["outdir"]
IMAGE_DIR = os.path.join(OUTDIR, "image.dir")
os.makedirs(OUTDIR, exist_ok=True)


# Build image-specific task information.
task_dict = {}

for slide, slide_config in config["slides"].items():
    for smp, fovlst in slide_config["samples"].items():
        assert smp not in task_dict, f"Duplicate sample name: {smp}."

        task_dict[smp] = {
            "slide": slide,
            "fovpos_csv": slide_config["fovpos_csv"],
            "morph_pfx": slide_config["morph_pfx"],
            "fovlst": fovlst,
        }

task_df = pd.DataFrame.from_dict(task_dict, orient="index")
task_df.to_csv(
    os.path.join(OUTDIR, "task.summary.csv"),
    index_label="sample",
)

SAMPLE_LST = sorted(task_dict)
SAMPLE_REGEX = "|".join(re.escape(smp) for smp in SAMPLE_LST)

wildcard_constraints:
    smp=f"(?:{SAMPLE_REGEX})"


# Implementation-specific image directory.
IMAGE_SOURCE = "ashlar.dir" if config["run_ashlar"] else "assembly.dir"
IMAGE_SOURCE_DIR = os.path.join(IMAGE_DIR, IMAGE_SOURCE)


rule full:
    input:
        expand(
            os.path.join(IMAGE_DIR, "{smp}", "FOV.positions.csv"),
            smp=SAMPLE_LST,
        ),
        expand(
            os.path.join(IMAGE_DIR, "{smp}", "image.ome.tiff"),
            smp=SAMPLE_LST,
        ),
    resources:
        **RESOURCES,
    params:
        lock=LOCK,
        imagedir=IMAGE_DIR,
    shell:
        """
        if [ {params.lock} = true ]; then
            chmod -R a-w {params.imagedir}
        fi

        echo "=======> All images done! <======="
        """


rule link_image:
    input:
        fov_csv=os.path.join(
            IMAGE_SOURCE_DIR,
            "{smp}",
            "FOV.positions.csv",
        ),
        morph_tiff=os.path.join(
            IMAGE_SOURCE_DIR,
            "{smp}",
            "image.ome.tiff",
        ),
    output:
        fov_csv=os.path.join(
            IMAGE_DIR,
            "{smp}",
            "FOV.positions.csv",
        ),
        morph_tiff=os.path.join(
            IMAGE_DIR,
            "{smp}",
            "image.ome.tiff",
        ),
    log:
        os.path.join(
            IMAGE_DIR,
            "{smp}",
            "cosmx_linkImage.log",
        ),
    resources:
        **RESOURCES,
    params:
        outdir=os.path.join(IMAGE_DIR, "{smp}"),
        source=IMAGE_SOURCE,
    shell:
        """
        mkdir -p {params.outdir}

        ln -sfn \
            ../{params.source}/{wildcards.smp}/FOV.positions.csv \
            {output.fov_csv}

        ln -sfn \
            ../{params.source}/{wildcards.smp}/image.ome.tiff \
            {output.morph_tiff}
        """


rule run_ashlar:
    input:
        os.path.join(
            IMAGE_DIR,
            "ashlar.dir",
            "{smp}",
            "grid_positions.csv",
        ),
    output:
        os.path.join(
            IMAGE_DIR,
            "ashlar.dir",
            "{smp}",
            "FOV.positions.csv",
        ),
        os.path.join(
            IMAGE_DIR,
            "ashlar.dir",
            "{smp}",
            "image.ome.tiff",
        ),
    log:
        os.path.join(
            IMAGE_DIR,
            "ashlar.dir",
            "{smp}",
            "cosmx_runAshlar.log",
        ),
    resources:
        **RESOURCES,
    params:
        outdir=os.path.join(
            IMAGE_DIR,
            "ashlar.dir",
            "{smp}",
        ),
        field_dir=os.path.join(
            IMAGE_DIR,
            "ashlar.dir",
            "{smp}",
            "field_links",
        ),
        name_pattern="F{series:05}.TIF",
        px_size=config["px_size"],
        overlap=config["ashlar"]["overlap"],
        max_shift=config["ashlar"]["max_shift"],
        ch_align=config["ashlar"]["align_channel"],
        ch_out=(
            f"--ch-out {config['channels_output']}"
            if config.get("channels_output")
            else ""
        ),
    shell:
        """
        python -m spatialhub.scripts.cosmx_runAshlar \
            --gridpos-csv {input} \
            --field-dir {params.field_dir} \
            --name-pattern '{params.name_pattern}' \
            --px-size {params.px_size} \
            --overlap {params.overlap} \
            --max-shift {params.max_shift} \
            --ch-align {params.ch_align} \
            {params.ch_out} \
            {params.outdir} >{log} 2>&1
        """


rule complete_grid:
    input:
        fov_pos=lambda wc: task_dict[wc.smp]["fovpos_csv"],
        mock_fov=os.path.join(
            IMAGE_DIR,
            "ashlar.dir",
            "mock_fov.tiff",
        ),
    output:
        directory(
            os.path.join(
                IMAGE_DIR,
                "ashlar.dir",
                "{smp}",
                "field_links",
            )
        ),
        os.path.join(
            IMAGE_DIR,
            "ashlar.dir",
            "{smp}",
            "grid_positions.csv",
        ),
    log:
        os.path.join(
            IMAGE_DIR,
            "ashlar.dir",
            "{smp}",
            "cosmx_completeGrid.log",
        ),
    resources:
        **RESOURCES,
    params:
        outdir=os.path.join(
            IMAGE_DIR,
            "ashlar.dir",
            "{smp}",
        ),
        fov_lst=lambda wc: task_dict[wc.smp]["fovlst"],
        m2d_pfx=lambda wc: task_dict[wc.smp]["morph_pfx"],
    shell:
        """
        python -m spatialhub.scripts.cosmx_completeGrid \
            --fov-csv {input.fov_pos} \
            --fov-lst {params.fov_lst} \
            --m2d-pfx {params.m2d_pfx} \
            --mock-tiff {input.mock_fov} \
            {params.outdir} >{log} 2>&1
        """


rule gen_blank_fov:
    output:
        os.path.join(
            IMAGE_DIR,
            "ashlar.dir",
            "mock_fov.tiff",
        ),
    log:
        os.path.join(
            IMAGE_DIR,
            "ashlar.dir",
            "cosmx_genBlankFOV.log",
        ),
    resources:
        **RESOURCES,
    params:
        width=config["fov_width"],
        height=config["fov_height"],
        n_channels=config["n_channels"],
    shell:
        """
        python -m spatialhub.scripts.cosmx_genBlankFOV \
            --width {params.width} \
            --height {params.height} \
            --n-channels {params.n_channels} \
            {output} >{log} 2>&1
        """


rule assemble_fovs:
    input:
        fov_pos=lambda wc: task_dict[wc.smp]["fovpos_csv"],
    output:
        os.path.join(
            IMAGE_DIR,
            "assembly.dir",
            "{smp}",
            "FOV.positions.csv",
        ),
        os.path.join(
            IMAGE_DIR,
            "assembly.dir",
            "{smp}",
            "image.ome.tiff",
        ),
    log:
        os.path.join(
            IMAGE_DIR,
            "assembly.dir",
            "{smp}",
            "cosmx_assembleFOVs.log",
        ),
    resources:
        **RESOURCES,
    params:
        outdir=os.path.join(
            IMAGE_DIR,
            "assembly.dir",
            "{smp}",
        ),
        fov_lst=lambda wc: task_dict[wc.smp]["fovlst"],
        morph_pfx=lambda wc: task_dict[wc.smp]["morph_pfx"],
        ch_out=(
            f"--ch-out {config['channels_output']}"
            if config.get("channels_output")
            else ""
        ),
    shell:
        """
        python -m spatialhub.scripts.cosmx_assembleFOVs \
            --fov-csv {input.fov_pos} \
            --fov-lst {params.fov_lst} \
            --m2d-pfx {params.morph_pfx} \
            {params.ch_out} \
            {params.outdir} >{log} 2>&1
        """