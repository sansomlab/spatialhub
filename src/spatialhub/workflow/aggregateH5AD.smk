import os


configfile: "aggregateH5AD.yaml"

# [NOTE] this key is set by the CLI, not the config file!
if "lock" in config:
    LOCK = config["lock"].lower()
else:
    LOCK = "false"

RESOURCES = {"threads": 4, "mem_mb": 16000, "time": "04:00:00", "partition": "short"}
RESOURCES.update(config.get("resources", {}))


# Extract relevant configuration and validate it
SAMPLES = sorted(z.removesuffix(".zarr") for z in os.listdir(config["zarr_dir"]) if z.endswith(".zarr"))
POINTS = [p.strip() for p in config["points_from"].split(",") if p.strip()]
SHAPES = [s.strip() for s in config["shapes_by"].split(",") if s.strip()]
COORDS = config.get("coords", "global")
AGGREG = config.get("agg_func", "sum")
CONCAT = config.get("concat", len(POINTS) > 1)

if not SAMPLES or not POINTS or not SHAPES:
    raise ValueError("Expected at least one sample, Points and Shapes element")
if len(set(POINTS)) != len(POINTS) or len(set(SHAPES)) != len(SHAPES):
    raise ValueError("Duplicate Points or Shapes elements specified in configuration")
if CONCAT and len(POINTS) <= 1:
    raise ValueError("Cannot run concat: please specify more than one Points element to concatenate")

# Define output H5AD file names for each Points table (trancripts, negative/control...), and for the combined table
def individual(wc):
    return os.path.join(config["h5ad_dir"], f"{wc.points}.{wc.shapes}.{COORDS}", f"{wc.cap}.h5ad")
def combined(wc):
    return os.path.join(config["h5ad_dir"], f"{'+'.join(POINTS)}.{wc.shapes}.{COORDS}", f"{wc.cap}.h5ad")


# Define individual and combined targets
INDIVIDUAL = expand(
    os.path.join(config["h5ad_dir"], "{points}.{shapes}." + COORDS, "{cap}.h5ad"),
    points=POINTS, shapes=SHAPES, cap=SAMPLES,
)
COMBINED = expand(
    os.path.join(config["h5ad_dir"], "+".join(POINTS) + ".{shapes}." + COORDS, "{cap}.h5ad"),
    shapes=SHAPES, cap=SAMPLES,
) if CONCAT and len(POINTS) > 1 else []

wildcard_constraints:
    points="|".join(POINTS)


rule full:
    input:
        INDIVIDUAL + COMBINED,
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
            "{points}.{shapes}." + COORDS, 
            "{cap}.h5ad",
        ),
    log:
        os.path.join(
            config["h5ad_dir"],
            "{points}.{shapes}." + COORDS, 
            "aggH5AD.{cap}.log",
        ),
    resources:
        **RESOURCES,
    params:
        coords=COORDS,
        agg_func=AGGREG,
    shell:
        """
        python -m spatialhub.scripts.aggregateH5AD \
            {output:q} \
            --in-zarr {input:q} \
            --points-from {wildcards.points:q} \
            --shapes-by {wildcards.shapes:q} \
            --coords {params.coords} \
            --agg-func {params.agg_func} \
            >{log} 2>&1
        """


rule concatenate_h5ads:
    input:
        lambda wc: [
            os.path.join(
                config["h5ad_dir"],
                f"{p}.{wc.shapes}.{COORDS}",
                f"{wc.cap}.h5ad"
                )
            for p in POINTS
        ]
    output:
        os.path.join(
            config["h5ad_dir"], 
            "+".join(POINTS) + ".{shapes}." + COORDS, 
            "{cap}.h5ad"
        )
    log:
        os.path.join(
            config["h5ad_dir"], 
            "+".join(POINTS) + ".{shapes}." + COORDS, 
            "concatH5AD.{cap}.log"
        )
    resources:
        **RESOURCES,
    params:
        points_from=config["points_from"],
        h5ad_dir=config["h5ad_dir"]
    shell:
        """
        python -m spatialhub.scripts.concatenateH5AD \
            {output:q} \
            {input:q} \
            >{log} 2>&1
        """