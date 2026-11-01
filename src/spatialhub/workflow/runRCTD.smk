import os

from spatialhub.workflow.utils import opt2cmd


configfile: "runRCTD.yaml"


# [NOTE] this key is set by the CLI, not the config file!
if "lock" in config:
    LOCK = config["lock"].lower()
else:
    LOCK = "false"

RESOURCES = {"threads": 4, "mem_mb": 16000, "time": "04:00:00", "partition": "short"}
RESOURCES.update(config.get("resources", {}))


task_dict = config["tasks"]


rule full:
    input:
        expand(os.path.join(config["outdir"], "{task}.rctd.pkl"), task=task_dict.keys()),
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


rule run_rctd:
    input:
        reference=lambda wc: task_dict[wc.task]["reference"],
        query=lambda wc: task_dict[wc.task]["query"],
    output:
        os.path.join(config["outdir"], "{task}.rctd.pkl"),
    log:
        os.path.join(config["outdir"], "{task}.rctd.log"),
    resources:
        **RESOURCES,
    params:
        celltype=lambda wc: task_dict[wc.task]["celltype"],
        mode=config.get("mode", "full"),
        xcoord_cmd=lambda wc: opt2cmd(task_dict[wc.task].get("xcoord"), "--xcoord"),
        ycoord_cmd=lambda wc: opt2cmd(task_dict[wc.task].get("ycoord"), "--ycoord"),
        min_umi_cmd=lambda wc: opt2cmd(task_dict[wc.task].get("min_umi"), "--min-umi"),
        min_umi_sigma_cmd=lambda wc: opt2cmd(
            task_dict[wc.task].get("min_umi_sigma"), "--min-umi-sigma"
        ),
    shell:
        """
        python -m spatialhub.scripts.runRCTD \
            {output} \
            --ref {input.reference} \
            --qry {input.query} \
            --celltype {params.celltype} \
            --mode {params.mode} \
            {params.xcoord_cmd} \
            {params.ycoord_cmd} \
            {params.min_umi_cmd} \
            {params.min_umi_sigma_cmd} \
            >{log} 2>&1
        """
