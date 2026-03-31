import os
import argparse
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import cnmf


if __name__ == "__main__":
    args = argparse.ArgumentParser(description="Prepare for cNMF parallel computation.")
    args.add_argument(
        "--anndata", required=True, help="Path to the input anndata object."
    )
    args.add_argument("--outdir", required=True, help="Path to the output directory.")
    args.add_argument(
        "--min_nsquares",
        default=100,
        help="Minimum number of squares expressing a gene.",
    )
    args.add_argument(
        "--hvgs",
        default=None,
        help="Path to the HVG table. By default, genes passing min_ncell filtration are used.",
    )
    args.add_argument(
        "--k_list",
        nargs="+",
        type=int,
        required=True,
        help="List of number of factors for cNMF.",
    )
    args.add_argument(
        "--n_iters", type=int, default=200, help="Number of NMF runs to perform [200]."
    )
    args.add_argument(
        "--X_corr", default=None, help="Path to the X_corr matrix [None]."
    )
    args = args.parse_args()

    outpath = os.path.join(args.outdir, "cnmf_tmp")
    if os.path.exists(outpath):
        raise FileExistsError(f"Temp folder {outpath} already exists.")

    adata = sc.read_h5ad(args.anndata)
    print("Original adata:\n", adata)
    sc.pp.filter_genes(adata, min_cells=args.min_nsquares)
    print(f"Filtered anndata (min_cells={args.min_nsquares}):\n", adata)
    if args.hvgs is not None:
        hvgs = pd.read_csv(args.hvgs, index_col=0)
        adata = adata[:, adata.var_names.isin(hvgs.index)].copy()
        print("Anndata with HVGs:\n", adata)

    if args.X_corr is not None:
        print("Loading X_corr from", args.X_corr)
        X_corr = np.load(args.X_corr)
        if X_corr.shape != adata.X.shape:
            raise ValueError(
                f"Shape of X_corr {X_corr.shape} does not match adata.X {adata.X.shape}."
            )
        adata.X = X_corr
    else:
        if "log1p" not in adata.layers:
            raise ValueError("adata.layers['log1p'] not found")
        if hasattr(adata.layers["log1p"], "data"):
            data = adata.layers["log1p"].data
        else:
            data = adata.layers["log1p"].ravel()
        if np.allclose(data, np.rint(data), atol=1e-6):
            raise ValueError("adata.layers['log1p'] looks like raw counts")
        adata.X = adata.layers["log1p"].toarray()

    print(np.min(adata.X), np.max(adata.X))

    if adata.shape[1] < max(args.k_list):
        raise ValueError(
            f"Number of genes ({adata.shape[1]}) is less than the maximum number of factors ({max(args.k_list)})."
        )
    if adata.shape[1] > 5000:
        warnings.warn(f"{adata.shape[1]} genes are used in factorisation.", UserWarning)
    adata.write_h5ad(os.path.join(args.outdir, "cnmf/adata.h5ad"))

    cnmf_obj = cnmf.cNMF(output_dir=args.outdir, name="cnmf")
    cnmf_obj.prepare(
        counts_fn=os.path.join(args.outdir, "cnmf/adata.h5ad"),
        components=args.k_list,
        num_highvar_genes=adata.shape[1],
        n_iter=args.n_iters,
        seed=0,
    )
