import os
import shlex


configfile: "addZarrShape.yaml"


# [NOTE] this key is set by the CLI, not the config file!
if "lock" in config:
    LOCK = config["lock"].lower()
else:
    LOCK = "false"

RESOURCES = {"threads": 4, "mem_mb": 16000, "time": "04:00:00", "partition": "short"}
RESOURCES.update(config.get("resources", {}))

SHAPE_DICT = config["shapes2add"]
ZARR_DIR = config["zarr_dir"]
SHAPE_KEY = config["shape_key"]
OUTFILE_SFX = os.path.join("shapes", SHAPE_KEY, "shapes.parquet")
GEOJSON_DIR = config["geojson_dir"]


rule full:
    input:
        expand(
            os.path.join(ZARR_DIR, "{smp}.zarr", OUTFILE_SFX),
            smp=SHAPE_DICT.keys(),
        ),
    resources:
        **RESOURCES,
    params:
        lock=str(LOCK).lower(),
    shell:
        """
        if [ "{params.lock}" = "true" ]; then
            chmod -R a-w {ZARR_DIR:q}/*.zarr
        fi
        echo "=======> All done! <======="
        """


rule prepare_geojson:
    input:
        lambda wc: SHAPE_DICT[wc.smp],
    output:
        os.path.join(GEOJSON_DIR, "{smp}", f"{SHAPE_KEY}.geojson"),
    log:
        os.path.join(GEOJSON_DIR, "{smp}", f"prep_geojson.{SHAPE_KEY}.log"),
    resources:
        **RESOURCES,
    params:
        t2g=lambda wc, input: input[0].lower().endswith((".tif", ".tiff")),
    shell:
        """
        if [ "{params.t2g}" = "True" ]; then
            python -m spatialhub.scripts.convTiff2Geojson \
                {output:q} \
                --in-tiff {input:q} \
                >{log:q} 2>&1
        else
            ln -s "$(realpath {input:q})" {output:q}
            echo "Input shape file is already a GeoJSON. Created symlink to {output}." >{log:q}
        fi
        """


rule add_shapes:
    input:
        shape_file=os.path.join(GEOJSON_DIR, "{smp}", f"{SHAPE_KEY}.geojson"),
        zarr_path=os.path.join(ZARR_DIR, "{smp}.zarr"),
    output:
        os.path.join(ZARR_DIR, "{smp}.zarr", OUTFILE_SFX),
    log:
        os.path.join(ZARR_DIR, f"add_shape.{{smp}}.{SHAPE_KEY}.log"),
    resources:
        **RESOURCES,
    params:
        shpdir=lambda wc, input: os.path.join(input.zarr_path, "shapes", SHAPE_KEY),
        shapekey=SHAPE_KEY,
        index=lambda wc: (
            "index"
            if SHAPE_DICT[wc.smp].lower().endswith((".tif", ".tiff"))
            else config.get("geojson_index", "index")
        ),
    shell:
        """
        rmdir {params.shpdir:q}
        python -m spatialhub.scripts.addZarrShape \
            {input.zarr_path:q} \
            --geojson-path {input.shape_file:q} \
            --key {params.shapekey:q} \
            --index {params.index:q} \
            >{log:q} 2>&1
        """
