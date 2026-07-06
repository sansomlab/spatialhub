import os
import pandas as pd
import warnings


configfile: "visiumhd_prepareBin2Cell.yaml"


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
            os.path.join(config["outdir"], "stardist/he.{mpp}-{pb_thres}.{fmt}"),
            mpp=map(float, config["mpp_lst"].split(",")),
            pb_thres=map(float, config["probability_threshold_lst"].split(",")),
            fmt=["tif", "npz"],
        ),
    resources:
        **RESOURCES,
    params:
        lock=LOCK,
        outdir=config["outdir"],
    shell:
        """
        if [ {params.lock} = true ]; then
            chmod -R a-w {params.outdir}
        fi
        echo "=======> All done! <======="
        """


rule bin2cell:
    input:
        sq02_dir=config["square_002um_dir"],
        src_img=config["source_image_path"],
        spa_dir=config["spaceranger_spatial_dir"],
    output:
        os.path.join(config["outdir"], "stardist/he.{mpp}-{pb_thres}.tif"),
        os.path.join(config["outdir"], "stardist/he.{mpp}-{pb_thres}.npz"),
    log:
        os.path.join(config["outdir"], "bin2cell.{mpp}-{pb_thres}.log"),
    resources:
        **RESOURCES,
    params:
        outdir=config["outdir"],
        mcells=config["min_cells"],
        mcounts=config["min_counts"],
        mpp="{mpp}",
        destripe_cmd=f"--destripe" if config.get("destripe", False) else "",
        pb_thres="{pb_thres}",
    shell:
        """
        python -m spatialhub.scripts.visiumhd_prepareBin2Cell \
            {params.outdir} \
            --sq2 {input.sq02_dir} \
            --sc-img {input.src_img} \
            --sr-spa {input.spa_dir} \
            --mcells {params.mcells} \
            --mcounts {params.mcounts} \
            --mpp {params.mpp} \
            {params.destripe_cmd} \
            --prob-thresh {params.pb_thres} \
            >{log} 2>&1
        """
