import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import ssam
import anndata

from skimage.exposure import rescale_intensity
from skimage.filters import gaussian, threshold_otsu
from skimage.morphology import remove_small_objects, remove_small_holes, closing, disk


def make_mask(
    img2pp,
    dna_weight=0.7,
    mem_weight=0.3,
    transform="off",
    min_obj_size=500,
    min_hole_size=1000,
    closing_dsk=3,
):
    # Combining DNA and Membrane channels to generate mask
    assert dna_weight + mem_weight == 1, "channel weights do not sum to 1."
    dna_ch = img2pp.sel(c="DNA").values
    dna_ch = rescale_intensity(dna_ch, in_range="image", out_range=(0, 1))
    mem_ch = img2pp.sel(c="Membrane").values
    mem_ch = rescale_intensity(mem_ch, in_range="image", out_range=(0, 1))
    fused = dna_weight * dna_ch + mem_weight * mem_ch

    # Non-linear transformation
    if transform == "log1p":
        fused = np.log1p(fused)
    elif transform == "sqrt":
        fused = np.sqrt(fused)
    elif transform != "off":
        assert False, f"unknown transformation method: {transform}"

    # Image filtering and compute threshold through Otsu's method
    fused_blur = gaussian(fused, sigma=1)
    thresh = threshold_otsu(fused_blur)

    # Making mask with morphological processing
    mask = fused_blur > thresh
    mask = remove_small_objects(mask, min_size=min_obj_size)
    mask = remove_small_holes(mask, area_threshold=min_hole_size)
    if closing_dsk > 0:
        mask = closing(mask, disk(closing_dsk))

    return fused, thresh, mask


def cmp_channels(img2pl):
    axes = plt.subplots(1, 3, figsize=(12, 5))[1]

    _, _, msk1 = make_mask(
        img2pl.copy(),
        dna_weight=1,
        mem_weight=0,
        transform="off",
        min_obj_size=0,
        min_hole_size=0,
        closing_dsk=0,
    )
    axes[0].imshow(msk1)
    axes[0].set_title("DNA only")

    _, _, msk2 = make_mask(
        img2pl.copy(),
        dna_weight=0,
        mem_weight=1,
        transform="off",
        min_obj_size=0,
        min_hole_size=0,
        closing_dsk=0,
    )
    axes[1].imshow(msk2)
    axes[1].set_title("Membrane only")

    _, _, msk3 = make_mask(
        img2pl.copy(),
        dna_weight=0.5,
        mem_weight=0.5,
        transform="off",
        min_obj_size=0,
        min_hole_size=0,
        closing_dsk=0,
    )
    axes[2].imshow(msk3)
    axes[2].set_title("0.5DNA + 0.5Membrane")

    plt.tight_layout()


def cmp_transfm(img2pl):
    axes = plt.subplots(2, 3, figsize=(12, 9), height_ratios=(2, 1))[1]

    fuse1, thres1, msk1 = make_mask(
        img2pl.copy(),
        dna_weight=0.5,
        mem_weight=0.5,
        transform="off",
        min_obj_size=0,
        min_hole_size=0,
        closing_dsk=0,
    )
    axes[0, 0].imshow(msk1)
    axes[0, 0].set_title("transformation off")
    sns.histplot(fuse1.flatten(), ax=axes[1, 0], bins=100).axvline(
        x=thres1, color="red", linestyle="dotted"
    )

    fuse2, thres2, msk2 = make_mask(
        img2pl.copy(),
        dna_weight=0.5,
        mem_weight=0.5,
        transform="log1p",
        min_obj_size=0,
        min_hole_size=0,
        closing_dsk=0,
    )
    axes[0, 1].imshow(msk2)
    axes[0, 1].set_title("log1p transformation")
    sns.histplot(fuse2.flatten(), ax=axes[1, 1], bins=100).axvline(
        x=thres2, color="red", linestyle="dotted"
    )

    fuse3, thres3, msk3 = make_mask(
        img2pl.copy(),
        dna_weight=0.5,
        mem_weight=0.5,
        transform="sqrt",
        min_obj_size=0,
        min_hole_size=0,
        closing_dsk=0,
    )
    axes[0, 2].imshow(msk3)
    axes[0, 2].set_title("sqrt transformation")
    sns.histplot(fuse3.flatten(), ax=axes[1, 2], bins=100).axvline(
        x=thres3, color="red", linestyle="dotted"
    )

    plt.tight_layout()


def cmp_technical(img2pl):
    axes = plt.subplots(2, 3, figsize=(12, 9))[1]

    _, _, msk11 = make_mask(
        img2pl.copy(),
        dna_weight=0.5,
        mem_weight=0.5,
        transform="sqrt",
        min_obj_size=200,
        min_hole_size=2000,
        closing_dsk=3,
    )
    axes[0, 0].imshow(msk11)
    axes[0, 0].set_title("obj>=200, hole>=2000, closing_disk=3")
    _, _, msk12 = make_mask(
        img2pl.copy(),
        dna_weight=0.5,
        mem_weight=0.5,
        transform="sqrt",
        min_obj_size=2000,
        min_hole_size=2000,
        closing_dsk=3,
    )
    axes[1, 0].imshow(msk12)
    axes[1, 0].set_title("obj>=2000, hole>=2000, closing_disk=3")

    _, _, msk21 = make_mask(
        img2pl.copy(),
        dna_weight=0.5,
        mem_weight=0.5,
        transform="sqrt",
        min_obj_size=1000,
        min_hole_size=500,
        closing_dsk=3,
    )
    axes[0, 1].imshow(msk21)
    axes[0, 1].set_title("obj>=1000, hole>=500, closing_disk=3")
    _, _, msk22 = make_mask(
        img2pl.copy(),
        dna_weight=0.5,
        mem_weight=0.5,
        transform="sqrt",
        min_obj_size=1000,
        min_hole_size=3000,
        closing_dsk=3,
    )
    axes[1, 1].imshow(msk22)
    axes[1, 1].set_title("obj>=1000, hole>=3000, closing_disk=3")

    _, _, msk31 = make_mask(
        img2pl.copy(),
        dna_weight=0.5,
        mem_weight=0.5,
        transform="sqrt",
        min_obj_size=1000,
        min_hole_size=2000,
        closing_dsk=1,
    )
    axes[0, 2].imshow(msk31)
    axes[0, 2].set_title("obj>=1000, hole>=2000, closing_disk=1")
    _, _, msk32 = make_mask(
        img2pl.copy(),
        dna_weight=0.5,
        mem_weight=0.5,
        transform="sqrt",
        min_obj_size=1000,
        min_hole_size=2000,
        closing_dsk=5,
    )
    axes[1, 2].imshow(msk32)
    axes[1, 2].set_title("obj>=1000, hole>=2000, closing_disk=5")

    plt.tight_layout()


def prepare_ssam_from_cosmx(dfpp, px2um=0.120280945):
    # Check existence of all required columns
    assert "x" in dfpp.columns, "dfpp does not have the required `x` column."
    assert "y" in dfpp.columns, "dfpp does not have the required `y` column."
    assert "gene" in dfpp.columns, "dfpp does not have the required `gene` column."

    # Backup original coordinates
    dfpp[["original_x_px", "original_y_px"]] = dfpp[["x", "y"]].copy()
    dfpp.x, dfpp.y = dfpp.x.astype(float), dfpp.y.astype(float)

    # Transform coordinates
    xmin, ymin, xmax, ymax = dfpp.x.min(), dfpp.y.min(), dfpp.x.max(), dfpp.y.max()
    dfpp.x -= xmin
    dfpp.y -= ymin
    dfpp.x *= px2um
    dfpp.y *= px2um

    # Get image width and height in micrometers
    um_width, um_height = dfpp.x.max(), dfpp.y.max()

    return (dfpp, xmin, ymin, xmax, ymax, um_width, um_height)


def analyse_ssam(
    dfpp,
    mdl_path=None,
    um_width=None,
    um_height=None,
    norm_thres=None,
    expr_thres=None,
    bandwidth=2.5,
    mask=None,
):
    assert mdl_path is not None, "Model storage path is required."
    assert isinstance(um_width, float), "um_width not correctly provided."
    assert isinstance(um_height, float), "um_height not correctly provided."

    # Initialise dataset and analysis objects
    ds = ssam.SSAMDataset(mdl_path)
    analysis = ssam.SSAMAnalysis(ds, verbose=True)

    # Run KDE and search local maxima
    analysis.run_kde(
        dfssam, width=um_width, height=um_height, bandwidth=bandwidth, re_run=True
    )
    if (norm_thres is not None) or (expr_thres is not None):
        analysis.set_thresholds(
            norm_threshold=norm_thres, expression_threshold=expr_thres
        )
    analysis.find_localmax(search_size=3, mask=mask)

    # Plot local maxima
    ds.plot_l1norm(cmap="Greys")
    ds.plot_localmax(c="Blue", s=0.2)
    plt.title(
        f"maxima by norm>{round(ds.norm_threshold, 3)}, expr>{round(ds.expression_threshold, 3)}"
    )

    # Normalise selected vectors and generate anndata
    analysis.normalize_vectors(normalize_vector=True, log_transform=True)
    adata = anndata.AnnData(
        ds.normalized_vectors,
        var=pd.DataFrame(index=ds.genes),
        obs=pd.DataFrame.from_dict({"x": ds.local_maxs[0], "y": ds.local_maxs[1]}),
    )

    return ds, analysis, adata
