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
        "--outdir", type=str, required=True, help="Path to save the RCTD results."
    )
    args = parser.parse_args()

    celltype_clean = re.sub(r"[^A-Za-z0-9_-]", "", args.celltype)
    opath = os.path.join(args.outdir, f"rctd_result.{args.mode}.{celltype_clean}.pkl")
    if os.path.exists(opath):
        raise FileExistsError(
            f"Output file {opath} already exists. Please rerun with a different path."
        )

    # Prepare data for RCTD
    ref = rctd.Reference(sc.read_h5ad(args.reference), cell_type_col=args.celltype)

    # Run RCTD
    cfg = rctd.RCTDConfig(UMI_min=1, UMI_min_sigma=1)
    result = rctd.run_rctd(sc.read_h5ad(args.query), ref, mode=args.mode, config=cfg)

    # Save results
    with open(opath, "wb") as f:
        pickle.dump(result, f)
