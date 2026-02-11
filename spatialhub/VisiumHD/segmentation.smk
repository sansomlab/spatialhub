import os


resources = config.get("resources", None)

spaceranger_dir = config["spaceranger_dir"]
fullres_image_dir = config["fullres_image_dir"]
samples = config["samples"]

outdir = config["output_dir"]

DEFAULT_PBTHRESH = [0.01]
DEFAULT_NMSTHRESH = [0.3]

mpp_values = config["mpp_values"]
best_mpp = config.get("best_mpp", None)
prob_thresh_values = config.get("prob_thresh_values", None)
if prob_thresh_values is None:
    prob_thresh_values = DEFAULT_PBTHRESH
nms_thresh_values = config.get("nms_thresh_values", None)
if nms_thresh_values is None:
    nms_thresh_values = DEFAULT_NMSTHRESH


rule full:
    input:
        expand(
            os.path.join(
                outdir, "{sample}", "stardist", "he.mpp{mpp}_pb{pb}_nms{nms}.npz"
            ),
            sample=samples,
            mpp=mpp_values,
            pb=DEFAULT_PBTHRESH,
            nms=DEFAULT_NMSTHRESH,
        ),
        expand(
            os.path.join(
                outdir, "{sample}", "stardist", "he.mpp{mpp}_pb{pb}_nms{nms}.npz"
            ),
            sample=samples,
            mpp=best_mpp,
            pb=prob_thresh_values,
            nms=nms_thresh_values,
        ),


rule preprocessing:
    input:
        fld=os.path.join(spaceranger_dir, "{sample}/outs/binned_outputs/square_002um"),
        img=os.path.join(fullres_image_dir, "{sample}.tif"),
    output:
        imgs=[
            os.path.join(outdir, "{sample}", "scaled_images", f"he.mpp{mpp}.tiff")
            for mpp in mpp_values
        ],
        h5ad=os.path.join(outdir, "{sample}", "preprocessed.h5ad"),
    log:
        os.path.join(outdir, "{sample}", "scale_he_image.{sample}.log"),
    params:
        smp="{sample}",
        minbins=config["preprocessing"].get("min_bins", 3),
        mincounts=config["preprocessing"].get("min_counts", 1),
        quantile=config["preprocessing"]["quantile"],
        mpp_values=config["preprocessing"]["mpp_values"],
    resources:
        **resources,
    run:
        import sys
        import scanpy as sc
        import bin2cell as b2c

        logfile = open(log[0], "w")
        sys.stdout = logfile
        sys.stderr = logfile

        adata = b2c.read_visium(os.path.join(input.fld), source_image_path=input.img)
        adata.var_names_make_unique()
        sc.pp.filter_genes(adata, min_cells=params.minbins)
        sc.pp.filter_cells(adata, min_counts=params.mincounts)
        b2c.destripe(adata, quantile=params.quantile)

        for mpp in params.mpp_values:
            spath = os.path.join(outdir, params.smp, f"scaled_images/he.mpp{mpp}.tiff")
            b2c.scaled_he_image(adata, mpp=mpp, save_path=spath)

        adata.write_h5ad(output.h5ad, compression="gzip")

        logfile.close()


rule segment_he:
    input:
        scaled_img=os.path.join(outdir, "{sample}", "scaled_images", "he.mpp{mpp}.tiff"),
    output:
        os.path.join(
            outdir, "{sample}", "stardist_he", "he.mpp{mpp}_pb{pb}_nms{nms}.npz"
        ),
    log:
        os.path.join(
            outdir, "{sample}", "stardist_he", "he.mpp{mpp}_pb{pb}_nms{nms}.log"
        ),
    params:
        sample="{sample}",
        mpp="{mpp}",
        prob_thresh="{pb}",
        nms_thresh="{nms}",
    resources:
        **resources,
    run:
        import gc
        import sys
        import bin2cell as b2c

        logfile = open(log[0], "w")
        sys.stdout = logfile
        sys.stderr = logfile

        img_dims = b2c.load_image(input.scaled_img).shape[:2]
        block_size = min(img_dims) / 2
        kwargs = {}
        if params.prob_thresh is not None:
            kwargs["prob_thresh"] = float(params.prob_thresh)
        if params.nms_thresh is not None:
            kwargs["nms_thresh"] = float(params.nms_thresh)
        gc.collect()

        b2c.stardist(
            image_path=input.scaled_img,
            block_size=block_size,
            labels_npz_path=output[0],
            stardist_model="2D_versatile_he",
            **kwargs,
        )

        logfile.close()


rule segment_gex:
    input:
        h5ad=os.path.join(outdir, "{sample}", "scaled_he.h5ad"),
    output:
        os.path.join(
            outdir, "{sample}", "stardist_gex", "gex.mpp{mpp}_pb{pb}_nms{nms}.npz"
        ),
    log:
        os.path.join(
            outdir, "{sample}", "stardist_gex", "gex.mpp{mpp}_pb{pb}_nms{nms}.log"
        ),
    params:
        sample="{sample}",
        mpp="{mpp}",
        prob_thresh="{pb}",
        nms_thresh="{nms}",
    resources:
        **resources,
    run:
        import gc
        import sys
        import bin2cell as b2c

        logfile = open(log[0], "w")
        sys.stdout = logfile
        sys.stderr = logfile

        adata = b2c.read_h5ad(input.h5ad)
        img_dims = adata.uns["image_shape"][:2]
        block_size = min(img_dims) / 2
        kwargs = {}
        if params.prob_thresh is not None:
            kwargs["prob_thresh"] = float(params.prob_thresh)
        if params.nms_thresh is not None:
            kwargs["nms_thresh"] = float(params.nms_thresh)
        gc.collect()

        b2c.stardist(
            adata,
            block_size=block_size,
            labels_npz_path=output[0],
            stardist_model="2D_versatile_he",
            **kwargs,
        )

        logfile.close()
