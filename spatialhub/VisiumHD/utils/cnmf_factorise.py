import argparse
import cnmf


if __name__ == "__main__":
    args = argparse.ArgumentParser(description="Prepare for cNMF parallel computation.")
    args.add_argument("--outdir", required=True, help="Path to the output directory.")
    args.add_argument(
        "--iworker", type=int, required=True, help="Worker ID for parallel execution."
    )
    args.add_argument(
        "--total_workers", type=int, required=True, help="Total number of workers."
    )
    args = args.parse_args()

    cnmf_obj = cnmf.cNMF(output_dir=args.outdir, name="cnmf")
    cnmf_obj.factorize(worker_i=args.iworker, total_workers=args.total_workers)
