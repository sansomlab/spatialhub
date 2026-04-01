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
            os.path.join(config["outdir"], "{tname}", "cnmf/{fname}"),
            tname=config["tasks"].keys(),
            fname=["cnmf.k_selection_stats.df.npz", "cnmf.k_selection.png"],
        ),
        # expand(
        #     os.path.join(OUTDIR, "{sample}", "cell2loc_nmf", "c2l_NMF.{nfactors}.pkl"),
        #     nfactors=config["nfact_list"],
        #     sample=SAMPLES,
        # ),


# rule cell2loc_nmf:
#     input:
#         anndata=os.path.join(ANNDATA_DIR, f"{{sample}}{ANNDATA_SUFFIX}"),
#     output:
#         os.path.join(OUTDIR, "{sample}", "cell2loc_nmf/c2l_NMF.{nfactors}.pkl"),
#     log:
#         os.path.join(OUTDIR, "{sample}", "cell2loc_nmf/c2l_NMF.{nfactors}.log"),
#     resources:
#         **RESOURCES,
#     params:
#         hvgs=os.path.join(HVGS_DIR, f"{{sample}}{HVGS_SUFFIX}"),
#         outdir=os.path.join(OUTDIR, "{sample}", "cell2loc_nmf"),
#         nfact="{nfactors}",
#         nrepeats=config["cell2location"]["n_repeats"],
#         sample_key=config["cell2location"].get("sample_key", None),
#         X_corr=config.get("X_corr", None),
#     shell:
#         """
#         if [ "{params.sample_key}" != "None" ]; then
#             params_sample_key="--sample_key {params.sample_key}"
#         else
#             params_sample_key=""
#         fi
#         if [ "{params.X_corr}" != "None" ]; then
#             params_X_corr="--X_corr {params.X_corr}"
#         else
#             params_X_corr=""
#         fi

#         python {UTILS_DIR}/cell2location_nmf.py \
#             --anndata {input.anndata} \
#             --hvgs {params.hvgs} \
#             --outdir {params.outdir} \
#             --nfactors {params.nfact} \
#             --nrepeats {params.nrepeats} \
#             $params_sample_key \
#             $params_X_corr \
#             > {log} 2>&1
#         """


rule cnmf_prepare:
    input:
        anndata=lambda wc: config["tasks"][wc.tname]["anndata"],
    output:
        os.path.join(config["outdir"], "{tname}", "cnmf/cnmf_tmp/cnmf.tpm.h5ad"),
        expand(
            os.path.join(config["outdir"], "{tname}", "cnmf/cnmf_tmp/{fname}"),
            tname=["{tname}"],
            fname=[
                "cnmf.tpm_stats.df.npz",
                "cnmf.norm_counts.h5ad",
                "cnmf.nmf_idvrun_params.yaml",
                "cnmf.nmf_params.df.npz",
            ],
        ),
    log:
        os.path.join(config["outdir"], "{tname}", "cnmf/cnmf_tmp/prepare.log"),
    resources:
        **RESOURCES,
    params:
        outdir=os.path.join(config["outdir"], "{tname}"),
        k_list=" ".join([str(k) for k in config["nfact_list"]]),
        min_nsquares=lambda wc: config["tasks"][wc.tname]["min_nsquares"],
        genes=lambda wc: config["tasks"][wc.tname].get("genes", None),
        n_iters=lambda wc: config["cnmf_niters"],
        X_corr=config.get("X_corr", None),
    shell:
        """
        if [ "{params.genes}" != "None" ]; then
            params_genes="--genes {params.genes}"
        else
            params_genes=""
        fi

        if [ "{params.X_corr}" != "None" ]; then
            params_X_corr="--X_corr {params.X_corr}"
        else
            params_X_corr=""
        fi

        python {UTILS_DIR}/cnmf_prepare.py \
            --anndata {input.anndata} \
            --min_nsquares {params.min_nsquares} \
            $params_genes \
            --outdir {params.outdir} \
            --k_list {params.k_list} \
            --n_iters {params.n_iters} \
            $params_X_corr \
            > {log} 2>&1
        """


rule cnmf_factorise:
    input:
        rules.cnmf_prepare.output,
    output:
        expand(
            os.path.join(
                config["outdir"],
                "{tname}",
                "cnmf/cnmf_tmp/cnmf.spectra.k_{k}.iter_{iw}.df.npz",
            ),
            tname=["{tname}"],
            k=config["nfact_list"],
            iw=["{iw}"],
        ),
    log:
        os.path.join(
            config["outdir"], "{tname}", "cnmf/cnmf_tmp/factorise.iw_start_{iw}.log"
        ),
    resources:
        **RESOURCES,
    params:
        outdir=os.path.join(config["outdir"], "{tname}"),
        iworker="{iw}",
        nworkers=config["cnmf_niters"],
    shell:
        """
        python {UTILS_DIR}/cnmf_factorise.py \
            --outdir {params.outdir} \
            --iworker {params.iworker} \
            --total_workers {params.nworkers} \
            > {log} 2>&1
        """


rule cnmf_combine:
    input:
        expand(
            rules.cnmf_factorise.output,
            iw=range(config["cnmf_niters"]),
            tname=["{tname}"],
        ),
    output:
        expand(
            os.path.join(config["outdir"], "{tname}", "cnmf/cnmf.{fname}"),
            tname=["{tname}"],
            fname=["k_selection_stats.df.npz", "k_selection.png"],
        ),
    resources:
        **RESOURCES,
    params:
        outdir=os.path.join(config["outdir"], "{tname}"),
    run:
        import cnmf

        cnmf_obj = cnmf.cNMF(outdir=params.outdir, name="cnmf")
        cnmf_obj.combine()
        cnmf_obj.k_selection_plot()
