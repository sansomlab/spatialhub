import os
import pandas as pd
import warnings

RESOURCES = {
    "threads": 1,
    "mem_mb": 4000,
    "time": "00:30:00",
    "ncpus": 1,
    "partition": "short",
}
RESOURCES.update(config.get("resources", {}))

UTILS_DIR = os.path.join(workflow.basedir, "utils")

metadata = pd.read_csv(config["meta_tsv"], sep="\t")
assert "sample_id" in metadata.columns, "sample_id column not found in metadata."
assert "slide_id" in metadata.columns, "slide_id column not found in metadata."
assert "fov_sequence" in metadata.columns, "fov_sequence column not found in metadata."
metadata.set_index("sample_id", inplace=True)

samples = metadata.index.unique()

slides = metadata["slide_id"].unique()
slide_dict = config["slide_dict"]
for slide in slides:
    assert slide in slide_dict, f"Slide {slide} not found in slide_dict."
    assert "flat_pfx" in slide_dict[slide], f"flat_pfx not found for slide {slide}."
    assert "morph_pfx" in slide_dict[slide], f"morph_pfx not found for slide {slide}."


rule full:
    input:
        expand(os.path.join(config["out_dir"], "{sid}.png"), sid=samples),
    resources:
        threads=1,
        mem_mb=1000,
        time="00:02:00",
        ncpus=1,
        partition=RESOURCES["partition"],
    params:
        out_dir=config["out_dir"],
    shell:
        """
        chmod -R a-w "{params.out_dir}"
        """


rule mkzarr_smp:
    output:
        outpath=directory(os.path.join(config["out_dir"], "{sid}.zarr")),
    log:
        os.path.join(config["out_dir"], "{sid}_makezarr.log"),
    resources:
        **RESOURCES,
    params:
        flat_pfx=lambda wc: os.path.join(
            slide_dict[metadata.loc[wc.sid, "slide_id"]]["flat_pfx"]
        ),
        include_polygon=config.get("include_polygon", False).__str__().lower(),
        ctrl_probe_regex=config.get("ctrl_probe_regex", "Negative|System|Blank"),
        morph_pfx=lambda wc: os.path.join(
            slide_dict[metadata.loc[wc.sid, "slide_id"]]["morph_pfx"]
        ),
        fov_lst=lambda wc: metadata.loc[wc.sid, "fov_sequence"],
        channel_lst=lambda wc: config.get("channel_lst", None),
        fov_width=config.get("fov_width", 4256),
        fov_height=config.get("fov_height", 4256),
        scale_factors=config.get("scale_factors", "2,4"),
    shell:
        """
        if [ "{params.include_polygon}" = "true" ]; then
            pg_csv_arg="--pg-csv {params.flat_pfx}-polygons.csv.gz"
        else
            pg_csv_arg=""
        fi

        if [ "{params.channel_lst}" == "None" ]; then
            channel_arg=""
        else
            channel_arg="--channel-lst {params.channel_lst}"
        fi

        python {UTILS_DIR}/make_zarr.py \
            --out-zarr "{output.outpath}" \
            --fovpos-csv "{params.flat_pfx}_fov_positions_file.csv.gz" \
            --tx-csv "{params.flat_pfx}_tx_file.csv.gz" \
            $pg_csv_arg \
            --ctrl-probe-regex "{params.ctrl_probe_regex}" \
            --morph-pfx "{params.morph_pfx}" \
            --fov-lst "{params.fov_lst}" \
            $channel_arg \
            --fov-width "{params.fov_width}" \
            --fov-height "{params.fov_height}" \
            --scale-factors "{params.scale_factors}" \
            > {log} 2>&1
        """


rule plot_check:
    input:
        os.path.join(config["out_dir"], "{sid}.zarr"),
    output:
        os.path.join(config["out_dir"], "{sid}.png"),
    log:
        os.path.join(config["out_dir"], "{sid}_checkimage.log"),
    resources:
        **RESOURCES,
    params:
        sid="{sid}",
    shell:
        """
        python {UTILS_DIR}/plot_zarr.py \
            --zarr "{input}" \
            --output "{output}" \
            > {log} 2>&1
        """
