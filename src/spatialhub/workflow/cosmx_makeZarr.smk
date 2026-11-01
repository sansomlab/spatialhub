import os
import pandas as pd


configfile: "cosmx_makeZarr.yaml"


# [NOTE] this key is set by the CLI, not the config file!
if "lock" in config:
    LOCK = config["lock"].lower()
else:
    LOCK = "false"

RESOURCES = {"threads": 4, "mem_mb": 16000, "time": "04:00:00", "partition": "short"}
RESOURCES.update(config.get("resources", {}))

OUTDIR = config["outdir"]
IMAGE_DIR = os.path.join(OUTDIR, "image.dir")
ZARR_DIR = os.path.join(OUTDIR, "zarr.dir")
os.makedirs(OUTDIR, exist_ok=True)

task_dict = {}
for slide, slide_config in config["slides"].items():
    tx_csv = slide_config["tx_csv"]
    pg_csv = slide_config.get("pg_csv", None)
    for smp, fovlst in slide_config["samples"].items():
        assert smp not in task_dict, f"Duplicate sample name: {smp}."
        task_dict[smp] = {
            "slide": slide,
            "tx_csv": tx_csv,
            "pg_csv": pg_csv,
            "fovlst": fovlst,
        }
task_df = pd.DataFrame.from_dict(task_dict, orient="index")
task_df.to_csv(os.path.join(config["outdir"], "task.summary.csv"), index_label="sample")

SAMPLE_LST = sorted(task_df.index.tolist())


rule full:
    input:
        expand(os.path.join(ZARR_DIR, "{smp}.zarr"), smp=SAMPLE_LST),
    resources:
        **RESOURCES,
    params:
        lock=LOCK,
        zarrdir=ZARR_DIR,
        imagedir=IMAGE_DIR,
    shell:
        """
        if [ {params.lock} = true ]; then
            chmod -R a-w {params.zarrdir} {params.imagedir}
        fi
        echo "=======> All done! <======="
        """


rule make_zarr:
    input:
        fov_csv=os.path.join(IMAGE_DIR, "{smp}", "FOV.positions.csv"),
        tx_csv=lambda wc: task_dict[wc.smp]["tx_csv"],
        morph_tiff=os.path.join(IMAGE_DIR, "{smp}", "image.ome.tiff"),
    output:
        directory(os.path.join(ZARR_DIR, "{smp}.zarr")),
    log:
        os.path.join(ZARR_DIR, "cosmx_makeZarr.{smp}.log"),
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
