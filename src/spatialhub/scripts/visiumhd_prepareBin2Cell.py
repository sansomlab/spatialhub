import os
import scanpy as sc
import numpy as np
import bin2cell as b2c
import matplotlib.pyplot as plt


from argparse import ArgumentParser as AP
from spatialhub.scripts.utils import die, print_arguments, RESET, GREEN


def main():
    p = AP(description="Perform Stardist steps for Bin2Cell Analysis.")
    p.add_argument("outdir", help="The output directory.")
    p.add_argument("--sq2", required=True, help="The SpaceRanger square_002um folder.")
    p.add_argument("--sc-img", required=True, help="The source image tif.")
    p.add_argument("--sr-spa", required=True, help="The SpaceRanger spatial folder.")
    p.add_argument("--mcells", type=int, default=3, help="Min square number per gene.")
    p.add_argument("--mcounts", type=int, default=1, help="Min counts per square.")
    p.add_argument("--mpp", type=float, default=0.5, help="Microns per pixel.")
    p.add_argument("--destripe", action="store_true", help="Destripe the image.")
    p.add_argument(
        "--prob-thresh",
        type=float,
        default=0.01,
        help="Stardist probability threshold.",
    )
    args = p.parse_args()
    print_arguments(args)

    out_tif = os.path.join(
        args.outdir, "stardist", f"he.{args.mpp}-{args.prob_thresh}.tif"
    )
    if os.path.exists(out_tif):
        raise FileExistsError(f"output files {out_tif} already exist")
    out_npz = os.path.join(
        args.outdir, "stardist", f"he.{args.mpp}-{args.prob_thresh}.npz"
    )
    if os.path.exists(out_npz):
        raise FileExistsError(f"output files {out_npz} already exist")
    out_dir = os.path.dirname(out_tif)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    adata = b2c.read_visium(
        args.sq2,
        source_image_path=args.sc_img,
        spaceranger_image_path=args.sr_spa,
    )
    adata.var_names_make_unique()

    sc.pp.filter_genes(adata, min_cells=args.mcells)
    sc.pp.filter_cells(adata, min_counts=args.mcounts)

    b2c.scaled_he_image(adata, mpp=args.mpp, save_path=out_tif)

    if args.destripe:
        b2c.destripe(adata)

    b2c.stardist(
        image_path=out_tif,
        labels_npz_path=out_npz,
        stardist_model="2D_versatile_he",
        prob_thresh=args.prob_thresh,
    )

    print(f"{GREEN}Successfully prepared for Bin2Cell Analysis.{RESET}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        die("Failed to create Zarr dataset.")
