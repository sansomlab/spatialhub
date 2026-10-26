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

# Points elements used to generate individual and concatenated H5ADs (if applicable)
POINTS_LIST = list(map(str.strip, config["points_from"].split(",")))
POINTS_LABEL = "+".join(POINTS_LIST)


if config["concat"]:
    FULL_TARGETS = expand(
        os.path.join(
            config["h5ad_dir"],
            f"{POINTS_LABEL}.{{shapes}}.{{agg}}.{{coords}}",
            "{cap}.h5ad",
        ),
        cap=CAPTURE_IDS,
        shapes=map(str.strip, config["shapes_by"].split(",")),
        agg=map(str.strip, config["agg_func"].split(",")),
        coords=map(str.strip, config["coords"].split(",")),
    )
else:
    FULL_TARGETS = expand(
        os.path.join(
            config["h5ad_dir"],
            "{points}.{shapes}.{agg}.{coords}",
            "{cap}.h5ad",
        ),
        cap=CAPTURE_IDS,
        points=map(str.strip, config["points_from"].split(",")),
        shapes=map(str.strip, config["shapes_by"].split(",")),
        agg=map(str.strip, config["agg_func"].split(",")),
        coords=map(str.strip, config["coords"].split(",")),
    )

POINTS = "|".join(map(str.strip, config["points_from"].split(",")))
wildcard_constraints:
    points=POINTS


rule full:
    input:
        FULL_TARGETS,
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
        os.path.join(
            config["h5ad_dir"],
            "{points}.{shapes}.{agg}.{coords}",
            "{cap}.h5ad",
        ),
    log:
        os.path.join(
            config["h5ad_dir"],
            "{points}.{shapes}.{agg}.{coords}",
            "aggH5AD.{cap}.log",
        ),
    resources:
        **RESOURCES,
    params:
        cap="{cap}",
        points="{points}",
        shapes="{shapes}",
        agg="{agg}",
        coords="{coords}",
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


rule concatenate_h5ads:
    input:
        lambda wc: expand(
            os.path.join(
                config["h5ad_dir"],
                "{points}.{shapes}.{agg}.{coords}",
                "{cap}.h5ad",
            ),
            cap=[wc.cap],
            points=POINTS_LIST,
            shapes=[wc.shapes],
            agg=[wc.agg],
            coords=[wc.coords],
        ),
    output:
        os.path.join(
            config["h5ad_dir"],
            f"{POINTS_LABEL}.{{shapes}}.{{agg}}.{{coords}}",
            "{cap}.h5ad",
        ),
    log:
        os.path.join(
            config["h5ad_dir"],
            f"{POINTS_LABEL}.{{shapes}}.{{agg}}.{{coords}}",
            "aggH5AD.{cap}.log",
        ),
    resources:
        **RESOURCES,
    params:
        points_from=config["points_from"],
        h5ad_dir=config["h5ad_dir"]
    shell:
        """
        python -m spatialhub.scripts.concatenateH5AD \
            {output} \
            --h5ad_dir "{params.h5ad_dir}" \
            --sample_id "{wildcards.cap}" \
            --points-from "{params.points_from}" \
            --shapes "{wildcards.shapes}" \
            --coords "{wildcards.coords}" \
            --agg-func "{wildcards.agg}" \
            >{log} 2>&1
        """