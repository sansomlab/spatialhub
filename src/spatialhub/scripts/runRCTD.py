import os
import rctd
import scanpy as sc
import pickle

from argparse import ArgumentParser as AP
from spatialhub.scripts.utils import die, GREEN, RESET, print_arguments

MODE_LST = ["full", "multi", "doublet"]


def main():
    p = AP(description="Deconvolute cell abundances with RCTD.")
    p.add_argument("outpath", help="Path to save the RCTD model.")
    p.add_argument("--ref", required=True, help="H5AD for the single-cell reference.")
    p.add_argument("--qry", required=True, help="H5AD for the spatial data in query.")
    p.add_argument("--celltype", required=True, help="Cell type column for reference.")
    p.add_argument("--mode", default="full", choices=MODE_LST, help="RCTD mode.")
    p.add_argument("--xcoord", default="x", help="X coordinates column in query data.")
    p.add_argument("--ycoord", default="y", help="Y coordinates column in query data.")
    p.add_argument("--min-umi", type=int, default=100, help="Minimum per pixel UMI.")
    p.add_argument(
        "--min-umi-sigma",
        type=int,
        default=300,
        help="Minimum per pixel UMI to estimate sigma.",
    )
    args = p.parse_args()
    print_arguments(args)

    if not args.outpath.endswith(".pkl"):
        raise ValueError("Output path must end with '.pkl'")
    if os.path.exists(args.outpath):
        raise FileExistsError(f"Output file {args.outpath} already exists")
    out_dir = os.path.dirname(args.outpath)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # Load reference data
    rdata = sc.read_h5ad(args.ref)
    print(rdata, "\n", rdata.X[rdata.X > 5], "\n")
    if "counts" in rdata.layers:
        rdata.X = rdata.layers["counts"].copy()
    else:
        raise ValueError("Reference data must contain a 'counts' layer.")

    if args.celltype not in rdata.obs.columns:
        raise ValueError(f"Column '{args.celltype}' not found in reference data")

    # Load query data
    qdata = sc.read_h5ad(args.qry)
    print(qdata, "\n", qdata.X[qdata.X > 5], "\n")
    if "counts" in qdata.layers:
        qdata.X = qdata.layers["counts"].copy()
    else:
        raise ValueError("Query data must contain a 'counts' layer.")

    missing = {args.xcoord, args.ycoord} - set(qdata.obs.columns)
    if missing:
        raise ValueError(f"Missing columns in query data: {missing}")

    if args.xcoord != "x":
        print(f"Copying {args.xcoord} column to 'x' for RCTD compatibility.")
        qdata.obs["x"] = qdata.obs[args.xcoord].values
    if args.ycoord != "y":
        print(f"Copying {args.ycoord} column to 'y' for RCTD compatibility.")
        qdata.obs["y"] = qdata.obs[args.ycoord].values

    # Create RCTD reference
    ref = rctd.Reference(rdata, cell_type_col=args.celltype)
    print(f"Reference profiles: {ref.profiles.shape}")
    print(f"Cell types: {ref.cell_type_names}")

    # Run RCTD
    cfg = rctd.RCTDConfig(UMI_min=args.min_umi, UMI_min_sigma=args.min_umi_sigma)
    result = rctd.run_rctd(qdata, ref, mode=args.mode, config=cfg)

    # Save results
    with open(args.outpath, "wb") as f:
        pickle.dump(result, f)

    print(f"{GREEN}RCTD completed successfully.{RESET}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        die("Failed to run RCTD.")
