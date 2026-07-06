import os

from spatialhub.workflow.utils import opt2cmd


configfile: "shared_runCell2Location.yaml"


# [NOTE] this key is set by the CLI, not the config file!
if "lock" in config:
    LOCK = config["lock"].lower()
else:
    LOCK = "false"

RESOURCES = {"threads": 4, "mem_mb": 16000, "time": "04:00:00", "partition": "short"}
RESOURCES.update(config.get("resources", {}))


rule full:
    input:
        os.path.join(config["outdir"], "mdl.reg.pt"),
        os.path.join(config["outdir"], "inf_aver.csv"),
        os.path.join(config["outdir"], "mdl.c2l.pt"),
        expand(
            os.path.join(config["outdir"], "{qt}_cell_abundance_w_sf.csv"),
            qt=map(str.strip, config["qt_lst"].split(",")),
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


QT_LST = map(str.strip, config["qt_lst"].split(","))


rule run_cell2location:
    input:
        reference=config["reference"],
        query=config["query"],
    output:
        os.path.join(config["outdir"], "mdl.reg.pt"),
        os.path.join(config["outdir"], "inf_aver.csv"),
        os.path.join(config["outdir"], "mdl.c2l.pt"),
        expand(
            os.path.join(config["outdir"], "{qt}_cell_abundance_w_sf.csv"), qt=QT_LST
        ),
    log:
        os.path.join(config["outdir"], "cell2location.log"),
    resources:
        **RESOURCES,
    params:
        outdir=config["outdir"],
        label_key=config["label_key"],
        batch_key_cmd=opt2cmd(config.get("batch_key"), "--batch-key"),
        cat_cmd=opt2cmd(config.get("categorical_covariates"), "--cat-keys"),
        ctn_cmd=opt2cmd(config.get("continuous_covariates"), "--ctn-keys"),
        mepoch_reg_cmd=opt2cmd(config.get("max_epochs_reg"), "--max-epochs-reg"),
        qt_lst_cmd=opt2cmd(config.get("qt_lst"), "--qt-lst"),
        n_cellperloc_cmd=opt2cmd(config.get("n_cellperloc"), "--n-cellperloc"),
        detect_alpha_cmd=opt2cmd(config.get("detection_alpha"), "--detection-alpha"),
        mepoch_c2l_cmd=opt2cmd(config.get("max_epochs_c2l"), "--max-epochs-c2l"),
    shell:
        """
        python -m spatialhub.scripts.shared_runCell2Location \
            {params.outdir} \
            --ref {input.reference} \
            --qry {input.query} \
            --label-key {params.label_key} \
            {params.batch_key_cmd} \
            {params.cat_cmd} \
            {params.ctn_cmd} \
            {params.mepoch_reg_cmd} \
            {params.qt_lst_cmd} \
            {params.n_cellperloc_cmd} \
            {params.detect_alpha_cmd} \
            {params.mepoch_c2l_cmd} \
            >{log} 2>&1
        """
