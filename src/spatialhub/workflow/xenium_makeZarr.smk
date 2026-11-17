import os
import pandas as pd


configfile: "xenium_makeZarr.yaml"


# [NOTE] this key is set by the CLI, not the config file!
if "lock" in config:
    LOCK = config["lock"].lower()
else:
    LOCK = "false"

RESOURCES = {
    "mem_mb": 16000,
    "time": "04:00:00",
    "partition": "short",
}
RESOURCES.update(config.get("resources", {}))
THREADS = RESOURCES.pop("threads", 1)


slide_dict = config["slides"]
slide_df = pd.DataFrame.from_dict(slide_dict, orient="index")

os.makedirs(config["outdir"], exist_ok=True)
slide_df.to_csv(os.path.join(config["outdir"], "slide.summary.csv"), index_label="sample")

SLIDE_LST = list(slide_dict.keys())
# TBC: How to generate one Zarr per sample for multi-sample capture areas - prob similar to CosMx


IMPORT_ADATA = "--import-adata" if config.get("import_adata", False) else ""


rule full:
    input:
        expand(os.path.join(config["outdir"], "zarr.dir", "{slide}.zarr"), slide=SLIDE_LST),
    params:
        lock=LOCK,
        zarr_dir=os.path.join(config["outdir"], "zarr.dir"),
    resources:
        **RESOURCES,
    shell:
        """
        if [ {params.lock} = true ]; then
            chmod -R a-w {params.zarr_dir}
        fi
        echo "=======> All done! <======="
        """


rule make_zarr:
    input:
        xoa_dir=lambda wc: slide_dict[wc.slide]["xoa_dir"],
    output:
        outdir=directory(os.path.join(config["outdir"], "zarr.dir", "{slide}.zarr")),
    log:
        os.path.join(config["outdir"], "zarr.dir", "xenium_makeZarr.{slide}.log"),
    threads:
        THREADS,
    resources:
        **RESOURCES,
    params:
        import_adata=IMPORT_ADATA,
    shell:
        """
        python -m spatialhub.scripts.xenium_makeZarr \
            {output.outdir} \
            --xoa-dir {input.xoa_dir} \
            {params.import_adata} \
            >{log} 2>&1
        """
