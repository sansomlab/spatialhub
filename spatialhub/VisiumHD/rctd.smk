import os

RESOURCES = {
    "threads": 1,
    "mem_mb": 4000,
    "time": "00:30:00",
    "ncpus": 1,
    "partition": "short",
}
RESOURCES.update(config.get("resources", {}))

UTILS_DIR = os.path.join(workflow.basedir, "utils")


rule full:
    input:
        expand(
            os.path.join(config["outdir"], "rctd_results.{tname}.pkl"),
            tname=config["tasks"].keys(),
        ),
    params:
        outdir=config["outdir"],
    resources:
        **RESOURCES,
    shell:
        """
        chmod -R a-w {params.outdir}
        """


rule run_rctd:
    input:
        query=lambda wc: config["tasks"][wc.tname]["query"],
        reference=lambda wc: config["tasks"][wc.tname]["reference"],
    output:
        os.path.join(config["outdir"], "rctd_results.{tname}.pkl"),
    log:
        os.path.join(config["outdir"], "rctd.{tname}.log"),
    resources:
        **RESOURCES,
    params:
        outpath=os.path.join(config["outdir"], "rctd_results.{tname}.pkl"),
        xcoord=lambda wc: config["tasks"][wc.tname]["xcoord"],
        ycoord=lambda wc: config["tasks"][wc.tname]["ycoord"],
        celltype=lambda wc: config["tasks"][wc.tname]["celltype"],
        umi_min=lambda wc: config["tasks"][wc.tname]["umi_min"],
        umi_min_sigma=lambda wc: config["tasks"][wc.tname]["umi_min_sigma"],
    shell:
        """
        python {UTILS_DIR}/run_rctd.py \
            --reference {input.reference} \
            --query {input.query} \
            --outpath {params.outpath} \
            --celltype {params.celltype} \
            --umi_min {params.umi_min} \
            --umi_min_sigma {params.umi_min_sigma} \
            --xcoord {params.xcoord} \
            --ycoord {params.ycoord} \
            > {log} 2>&1
        """
