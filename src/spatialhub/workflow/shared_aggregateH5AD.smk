import os


configfile: "shared_aggregateH5AD.yaml"


# [NOTE] this key is set by the CLI, not the config file!
if "lock" in config:
    LOCK = config["lock"].lower()
else:
    LOCK = "false"

RESOURCES = {"threads": 4, "mem_mb": 16000, "time": "04:00:00", "partition": "short"}
RESOURCES.update(config.get("resources", {}))


rule full:
    input:
        expand(
            os.path.join(
                config["workdir"], "h5ad.dir", "{cap}.{points}.{shapes}.{agg}.h5ad"
            ),
            cap=map(str.strip, config["capture_ids"].split(",")),
            points=[config["points_from"]],
            shapes=[config["shapes_by"]],
            agg=[config["agg_func"]],
        ),
    resources:
        **RESOURCES,
    params:
        lock=LOCK,
        h5ad_dir=os.path.join(config["workdir"], "h5ad.dir"),
    shell:
        """
        if [ {params.lock} = true ]; then
            chmod -R a-w {params.h5ad_dir}
        fi
        echo "=======> All done! <======="
        """


rule make_h5ad:
    input:
        os.path.join(config["workdir"], "zarr.dir", "{cap}.zarr"),
    output:
        os.path.join(
            config["workdir"], "h5ad.dir", "{cap}.{points}.{shapes}.{agg}.h5ad"
        ),
    log:
        os.path.join(
            config["workdir"], "h5ad.dir", "{cap}.{points}.{shapes}.{agg}.aggH5AD.log"
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
        python -m spatialhub.scripts.shared_aggregateH5AD \
            {output} \
            --in-zarr {input} \
            --points-from {params.points} \
            --shapes-by {params.shapes} \
            --agg-func {params.agg} \
            --coords {params.coords} \
            >{log} 2>&1
        """
