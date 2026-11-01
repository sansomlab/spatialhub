import os


configfile: "aggregateH5AD.yaml"


# [NOTE] this key is set by the CLI, not the config file!
if "lock" in config:
    LOCK = config["lock"].lower()
else:
    LOCK = "false"

RESOURCES = {"threads": 4, "mem_mb": 16000, "time": "04:00:00", "partition": "short"}
RESOURCES.update(config.get("resources", {}))

CAPTURE_IDS = [
    cap.replace(".zarr", "")
    for cap in os.listdir(config["zarr_dir"])
    if cap.endswith(".zarr")
]


rule full:
    input:
        expand(
            os.path.join(
                config["h5ad_dir"],
                "{cap}.{points}.{shapes}.{agg}.{coords}.h5ad",
            ),
            cap=CAPTURE_IDS,
            points=map(str.strip, config["points_from"].split(",")),
            shapes=map(str.strip, config["shapes_by"].split(",")),
            agg=map(str.strip, config["agg_func"].split(",")),
            coords=map(str.strip, config["coords"].split(",")),
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
        os.path.join(config["h5ad_dir"], "{cap}.{points}.{shapes}.{agg}.{coords}.h5ad"),
    log:
        os.path.join(
            config["h5ad_dir"], "aggH5AD.{cap}.{points}.{shapes}.{agg}.{coords}.log"
        ),
    resources:
        **RESOURCES,
    params:
        cap="{cap}",
        points="{points}",
        shapes="{shapes}",
        agg="{agg}",
        coords=config["coords"],
    shell:
        """
        python -m spatialhub.scripts.aggregateH5AD \
            {output} \
            --in-zarr {input} \
            --points-from {params.points} \
            --shapes-by {params.shapes} \
            --agg-func {params.agg} \
            --coords {params.coords} \
            >{log} 2>&1
        """
