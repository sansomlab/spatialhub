import os
import pandas as pd


configfile: "cosmx_makeZarr.yaml"


LOCK = config["lock"].lower()  # [NOTE] this key is set by the CLI, not the config file!

RESOURCES = {"threads": 4, "mem_mb": 16000, "time": "04:00:00", "partition": "short"}
RESOURCES.update(config.get("resources", {}))

task_dict = {}
for slide, slide_config in config["slides"].items():
    fovpos_csv = slide_config["fovpos_csv"]
    tx_csv = slide_config["tx_csv"]
    morph_pfx = slide_config["morph_pfx"]
    pg_csv = slide_config.get("pg_csv", None)
    for smp, fovlst in slide_config["samples"].items():
        assert smp not in task_dict, f"Duplicate sample name: {smp}."
        task_dict[smp] = {
            "slide": slide,
            "fovpos_csv": fovpos_csv,
            "tx_csv": tx_csv,
            "morph_pfx": morph_pfx,
            "pg_csv": pg_csv,
            "fovlst": fovlst,
        }
task_df = pd.DataFrame.from_dict(task_dict, orient="index")
task_df.to_csv(os.path.join(config["outdir"], "task.summary.csv"), index_label="sample")

SAMPLE_LST = sorted(task_df.index.tolist())


if config["run_ashlar"]:
    zarr_indir = os.path.join(config["outdir"], "ashlar.dir")
else:
    zarr_indir = os.path.join(config["outdir"], "assembly.dir")


rule full:
    input:
        expand(os.path.join(config["outdir"], "zarr.dir", "{smp}.zarr"), smp=SAMPLE_LST),
    resources:
        **RESOURCES,
    params:
        lock=LOCK,
        zarrdir=os.path.join(config["outdir"], "zarr.dir"),
        zarr_indir=zarr_indir,
    shell:
        """
        if [ {params.lock} = true ]; then
            chmod -R a-w {params.zarrdir} {params.zarr_indir}
        fi
        echo "=======> All done! <======="
        """


rule make_zarr:
    input:
        fov_csv=os.path.join(zarr_indir, "{smp}", "FOV.positions.csv"),
        tx_csv=lambda wc: task_dict[wc.smp]["tx_csv"],
        morph_tiff=os.path.join(zarr_indir, "{smp}", "image.ome.tiff"),
    output:
        directory(os.path.join(config["outdir"], "zarr.dir", "{smp}.zarr")),
    log:
        os.path.join(config["outdir"], "zarr.dir", "cosmx_makeZarr.{smp}.log"),
    resources:
        **RESOURCES,
    params:
        pg_csv=lambda wc: (
            f"--pg-csv {task_dict[wc.smp]['pg_csv']}"
            if task_dict[wc.smp].get("pg_csv")
            else ""
        ),
        chnames=(
            f"--channels {config['channel_names']}"
            if config.get("channel_names")
            else ""
        ),
        ctrl_regex=(
            f"--ctrl-regex '{config['ctrl_regex']}'"
            if config.get("ctrl_regex")
            else ""
        ),
        other_regex=(
            f"--other-regex '{config['other_regex']}'"
            if config.get("other_regex")
            else ""
        ),
        scales=(f"--scales {config['scales']}" if config.get("scales") else ""),
    shell:
        """
        python -m spatialhub.scripts.cosmx_makeZarr \
            --fov-csv {input.fov_csv} \
            --tx-csv {input.tx_csv} \
            --morph-tiff {input.morph_tiff} \
            {params.pg_csv} \
            {params.chnames} \
            {params.ctrl_regex} \
            {params.other_regex} \
            {params.scales} \
            {output} >{log} 2>&1
        """


rule run_ashlar:
    input:
        os.path.join(config["outdir"], "ashlar.dir", "{smp}", "grid_positions.csv"),
    output:
        os.path.join(config["outdir"], "ashlar.dir", "{smp}", "FOV.positions.csv"),
        os.path.join(config["outdir"], "ashlar.dir", "{smp}", "image.ome.tiff"),
    log:
        os.path.join(config["outdir"], "ashlar.dir", "{smp}", "cosmx_runAshlar.log"),
    resources:
        **RESOURCES,
    params:
        outdir=os.path.join(config["outdir"], "ashlar.dir", "{smp}"),
        field_dir=os.path.join(config["outdir"], "ashlar.dir", "{smp}", "field_links"),
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
        mock_fov=os.path.join(config["outdir"], "ashlar.dir", "mock_fov.tiff"),
    output:
        directory(os.path.join(config["outdir"], "ashlar.dir", "{smp}", "field_links")),
        os.path.join(config["outdir"], "ashlar.dir", "{smp}", "grid_positions.csv"),
    log:
        os.path.join(config["outdir"], "ashlar.dir", "{smp}", "cosmx_completeGrid.log"),
    resources:
        **RESOURCES,
    params:
        outdir=os.path.join(config["outdir"], "ashlar.dir", "{smp}"),
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
        os.path.join(config["outdir"], "ashlar.dir", "mock_fov.tiff"),
    log:
        os.path.join(config["outdir"], "ashlar.dir", "cosmx_genBlankFOV.log"),
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
        os.path.join(config["outdir"], "assembly.dir", "{smp}", "FOV.positions.csv"),
        os.path.join(config["outdir"], "assembly.dir", "{smp}", "image.ome.tiff"),
    log:
        os.path.join(
            config["outdir"], "assembly.dir", "{smp}", "cosmx_assembleFOVs.log"
        ),
    resources:
        **RESOURCES,
    params:
        outdir=os.path.join(config["outdir"], "assembly.dir", "{smp}"),
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
