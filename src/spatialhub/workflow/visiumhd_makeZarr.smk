import os
import pandas as pd
import warnings


configfile: "visiumhd_makeZarr.yaml"


LOCK = config["lock"].lower()  # [NOTE] this key is set by the CLI, not the config file!

RESOURCES = {"threads": 4, "mem_mb": 16000, "time": "04:00:00", "partition": "short"}
RESOURCES.update(config.get("resources", {}))

task_dict = config["tasks"]
task_df = pd.DataFrame.from_dict(task_dict, orient="index")
task_df.to_csv(os.path.join(config["outdir"], "task.summary.csv"), index_label="sample")

CAP_LST = list(task_dict.keys())

binsizes = [2, 8, 16]
custom_bin_size = config["spaceranger"].get("custom_bin_size")
if custom_bin_size:
    binsizes += [custom_bin_size]
    custbsize_cmd = f"--custom-bin-size {custom_bin_size}"
else:
    custbsize_cmd = ""
fmtbinsizes = [f"{bsize:03}" for bsize in binsizes]

probeset = config["spaceranger"].get("probeset")
if probeset:
    probeset_cmd = f"--probe-set {probeset}"
else:
    probeset_cmd = ""


rule full:
    input:
        expand(os.path.join(config["outdir"], "zarr.dir", "{cap}.zarr"), cap=CAP_LST),
    resources:
        **RESOURCES,
    params:
        lock=LOCK,
        sranger_dir=os.path.join(config["outdir"], "spaceranger.dir"),
        zarr_dir=os.path.join(config["outdir"], "zarr.dir"),
    shell:
        """
        if [ {params.lock} = true ]; then
            chmod -R a-w {params.sranger_dir} {params.zarr_dir}
        fi
        echo "=======> All done! <======="
        """


rule make_zarr:
    input:
        os.path.join(config["outdir"], "spaceranger.dir", "{cap}"),
    output:
        directory(os.path.join(config["outdir"], "zarr.dir", "{cap}.zarr")),
    log:
        os.path.join(config["outdir"], "zarr.dir", "visiumhd_makeZarr.{cap}.log"),
    resources:
        **RESOURCES,
    params:
        outdir=os.path.join(config["outdir"], "zarr.dir"),
        cap="{cap}",
        use_filtered=config["make_zarr"].get("use_filtered", True),
    run:
        import warnings
        import spatialdata_io as spdio
        from contextlib import redirect_stdout, redirect_stderr

        with open(log[0], "w") as log_file:
            with redirect_stdout(log_file), redirect_stderr(log_file):
                with warnings.catch_warnings():
                    warnings.filterwarnings(
                        "ignore",
                        category=UserWarning,
                        message=".*Variable names are not unique.*",
                    )
                    sdata = spdio.visium_hd(
                        path=os.path.join(os.path.dirname(input[0]), "outs"),
                        dataset_id=params.cap,
                        filtered_counts_file=params.use_filtered,
                        load_segmentations_only=False,
                        load_nucleus_segmentations=True,
                        bins_as_squares=True,
                        annotate_table_by_labels=True,
                        fullres_image_file=task_dict[params.cap].get("Image"),
                        load_all_images=True,
                        var_names_make_unique=True,
                    )
                sdata.write(os.path.join(params.outdir, f"{params.cap}.zarr"))


rule spaceranger_count:
    output:
        outdir=directory(os.path.join(config["outdir"], "spaceranger.dir", "{cap}")),
    log:
        os.path.join(config["outdir"], "spaceranger.dir", "spaceranger.{cap}.log"),
    threads: RESOURCES["threads"]
    resources:
        **RESOURCES,
    params:
        transcriptome=config["spaceranger"]["transcriptome"],
        probeset_cmd=probeset_cmd,
        slidefiledir=os.path.join(config["spaceranger"]["slidefiledir"], ""),
        createbam=str(config["spaceranger"].get("createbam", False)).lower(),
        custbsize_cmd=custbsize_cmd,
        cap="{cap}",
        fastqdir=lambda wc: task_dict[wc.cap]["fq_dir"],
        cytaimage=lambda wc: task_dict[wc.cap]["cyta_image"],
        image_cmd=lambda wc: (
            f"--image {task_dict[wc.cap].get('ms_image')}"
            if task_dict[wc.cap].get("ms_image")
            else ""
        ),
        slide=lambda wc: task_dict[wc.cap]["slide_serial"],
        area=lambda wc: task_dict[wc.cap]["capture_area"],
    shell:
        """
        module load SpaceRanger

        spaceranger count \
            --id "{params.cap}" \
            --transcriptome "{params.transcriptome}" \
            {params.probeset_cmd} \
            --create-bam "{params.createbam}" \
            --output-dir {output.outdir} \
            --fastqs "{params.fastqdir}" \
            --cytaimage "{params.cytaimage}" \
            {params.image_cmd} \
            --slide "{params.slide}" \
            --slidefile "{params.slidefiledir}{params.slide}.vlf" \
            --area "{params.area}" \
            {params.custbsize_cmd} \
            >{log} 2>&1
        """
