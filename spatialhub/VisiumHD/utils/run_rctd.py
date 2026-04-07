import os
import re
import rctd
import argparse
import scanpy as sc
import pickle


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run RCTD for spatial deconvolution.")
    parser.add_argument(
        "--reference",
        required=True,
        help="Path to the reference single-cell data (h5ad file).",
    )
    parser.add_argument(
        "--query",
        type=str,
        required=True,
        help="Path to the spatial data for query (h5ad file).",
    )
    parser.add_argument(
        "--celltype", required=True, help="Column in reference data for cell types."
    )
    parser.add_argument(
        "--outpath", type=str, required=True, help="Path to save the RCTD results."
    )
    parser.add_argument(
        "--mode", default="multi", help="RCTD mode: 'doublet', 'multi', or 'full'."
    )
    parser.add_argument(
        "--xcoord", default="x", help="Column name for x coordinates in query data."
    )
    parser.add_argument(
        "--ycoord", default="y", help="Column name for y coordinates in query data."
    )
    parser.add_argument(
        "--umi_min", type=int, default=100, help="Minimum UMI count per pixel."
    )
    parser.add_argument(
        "--umi_min_sigma",
        type=int,
        default=300,
        help="Minimum UMI for sigma estimation.",
    )
    args = parser.parse_args()

    celltype_clean = re.sub(r"[^A-Za-z0-9_-]", "", args.celltype)
    if not args.outpath.endswith(".pkl"):
        raise ValueError("Output path must end with '.pkl'")
    if os.path.exists(args.outpath):
        raise FileExistsError(
            f"Output file {args.outpath} already exists. Please rerun with a different path."
        )

    print("Arguments:")
    for arg in vars(args):
        print(f"\t- {arg}: {getattr(args, arg)}")

    # Load reference data
    rdata = sc.read_h5ad(args.reference)
    print(rdata, "\n", rdata.X[rdata.X > 5], "\n")

    # Load query data
    qdata = sc.read_h5ad(args.query)
    if "counts" in qdata.layers:
        qdata.X = qdata.layers["counts"].copy()
    if args.xcoord != "x":
        print(f"Copying {args.xcoord} column to 'x' for RCTD compatibility.")
        qdata.obs["x"] = qdata.obs[args.xcoord].values
    if args.ycoord != "y":
        print(f"Copying {args.ycoord} column to 'y' for RCTD compatibility.")
        qdata.obs["y"] = qdata.obs[args.ycoord].values
    print(qdata, "\n", qdata.X[qdata.X > 5], "\n")

    # Create RCTD reference
    ref = rctd.Reference(rdata, cell_type_col=args.celltype)
    print(f"Reference profiles: {ref.profiles.shape}")
    print(f"Cell types: {ref.cell_type_names}")

    # Run RCTD
    cfg = rctd.RCTDConfig(UMI_min=args.umi_min, UMI_min_sigma=args.umi_min_sigma)
    result = rctd.run_rctd(qdata, ref, mode=args.mode, config=cfg)

    # Save results
    with open(args.outpath, "wb") as f:
        pickle.dump(result, f)
