import os


configfile: "aggregateH5AD.yaml"


# [NOTE] this key is set by the CLI, not the config file!
if "lock" in config:
    LOCK = config["lock"].lower()
else:
    LOCK = "false"

RESOURCES = {"threads": 4, "mem_mb": 16000, "time": "04:00:00", "partition": "short"}
RESOURCES.update(config.get("resources", {}))

CAPTURE_IDS = sorted(
    cap.removesuffix(".zarr")
    for cap in os.listdir(config["zarr_dir"])
    if cap.endswith(".zarr")
)

BYSAMPLE_OUTPUT_DICT = config["bysample_outputs"]


rule full:
    input:
        expand(
            os.path.join(config["h5ad_dir"], "{cap}.{outname}.h5ad"),
            cap=CAPTURE_IDS,
            outname=list(BYSAMPLE_OUTPUT_DICT.keys()),
        ),
    resources:
        **RESOURCES,
    params:
        lock=LOCK,
        h5ad_dir=config["h5ad_dir"],
    shell:
        """
        if [ {params.lock} = true ]; then
            chmod -R a-w {params.h5ad_dir}
        fi
        echo "=======> All done! <======="
        """


rule make_h5ad:
    input:
        os.path.join(config["zarr_dir"], "{cap}.zarr"),
    output:
        os.path.join(config["h5ad_dir"], "{cap}.{outname}.h5ad"),
    log:
        os.path.join(config["h5ad_dir"], "aggH5AD.{cap}.{outname}.log"),
    resources:
        **RESOURCES,
    params:
        points=lambda wc: BYSAMPLE_OUTPUT_DICT[wc.outname]["points_from"],
        shapes=lambda wc: BYSAMPLE_OUTPUT_DICT[wc.outname]["shapes_by"],
        seg_idx=lambda wc: BYSAMPLE_OUTPUT_DICT[wc.outname].get("seg_idx"),
        agg_func=lambda wc: BYSAMPLE_OUTPUT_DICT[wc.outname].get("agg_func", "sum"),
        agg_name=lambda wc: BYSAMPLE_OUTPUT_DICT[wc.outname].get("agg_name", "counts"),
        coords=lambda wc: BYSAMPLE_OUTPUT_DICT[wc.outname].get("coords", "global"),
    shell:
        """
        if [ "{params.seg_idx}" = "None" ]; then
            seg_idx_arg=""
        else
            seg_idx_arg="--seg-idx {params.seg_idx}"
        fi

        python -m spatialhub.scripts.aggregateH5AD \
            {output} \
            --in-zarr {input} \
            --points-from {params.points} \
            --shapes-by {params.shapes} \
            $seg_idx_arg \
            --agg-func {params.agg_func} \
            --agg-name {params.agg_name} \
            --coords {params.coords} \
            >{log} 2>&1
        """
