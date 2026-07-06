import os


configfile: "shared_extractH5AD.yaml"


# [NOTE] this key is set by the CLI, not the config file!
if "lock" in config:
    LOCK = config["lock"].lower()
else:
    LOCK = "false"

RESOURCES = {"threads": 4, "mem_mb": 16000, "time": "04:00:00", "partition": "short"}
RESOURCES.update(config.get("resources", {}))

extract_tables = map(str.strip, config["extract_tables"].split(","))


rule full:
    input:
        expand(
            os.path.join(config["workdir"], "h5ad.dir", "{cap}.{table}.h5ad"),
            table=extract_tables,
            cap=map(str.strip, config["capture_ids"].split(",")),
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
        expand(
            os.path.join(config["workdir"], "h5ad.dir", "{cap}.{table}.h5ad"),
            table="{table}",
            cap="{cap}",
        ),
    log:
        os.path.join(config["workdir"], "h5ad.dir", "extractH5AD.{cap}.{table}.log"),
    resources:
        **RESOURCES,
    params:
        outdir=os.path.join(config["workdir"], "h5ad.dir"),
        cap="{cap}",
        tables="{table}",
    shell:
        """
        python -m spatialhub.scripts.shared_extractH5AD \
            {params.outdir}/{params.cap}.{params.tables}.h5ad \
            --in-zarr {input} \
            --extract {params.tables} \
            >{log} 2>&1
        """
