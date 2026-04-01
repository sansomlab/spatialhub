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
        "--celltype",
        required=True,
        help="Column in reference data for cell types.",
    )
    parser.add_argument(
        "--query",
        type=str,
        required=True,
        help="Path to the spatial data for query (h5ad file).",
    )
    parser.add_argument(
        "--mode",
        type=str,
        default="multi",
        help="RCTD mode: 'doublet', 'multi', or 'full'.",
    )
    parser.add_argument(
        "--xcoord", default="x", help="Column name for x coordinates in query data."
    )
    parser.add_argument(
        "--ycoord", default="y", help="Column name for y coordinates in query data."
    )
    parser.add_argument(
        "--outdir", type=str, required=True, help="Path to save the RCTD results."
    )
    args = parser.parse_args()

    celltype_clean = re.sub(r"[^A-Za-z0-9_-]", "", args.celltype)
    opath = os.path.join(args.outdir, f"rctd_result.{args.mode}.{celltype_clean}.pkl")
    if os.path.exists(opath):
        raise FileExistsError(
            f"Output file {opath} already exists. Please rerun with a different path."
        )

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
    cfg = rctd.RCTDConfig()
    result = rctd.run_rctd(qdata, ref, mode=args.mode, config=cfg)

    # Save results
    with open(opath, "wb") as f:
        pickle.dump(result, f)
