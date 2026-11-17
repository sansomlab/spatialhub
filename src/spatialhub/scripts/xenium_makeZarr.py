import os
import spatialdata_io as spdio

from argparse import ArgumentParser as AP
from spatialhub.scripts.utils import die, print_arguments, RESET, GREEN


def main():
    p = AP(description="Create Zarr dataset from Xenium data.")
    p.add_argument("output", help="Output Zarr directory.")
    p.add_argument("--xoa-dir", required=True, help="Path to the Xenium Onboard Analysis (XOA) pipeline directory.")
    p.add_argument("--import-adata", action="store_true", help="Import gene x cell data computed by Xenium Onboard Analysis software.")
    args = p.parse_args()
    print_arguments(args)

    if os.path.exists(args.output):
        raise FileExistsError(f"output file {args.output} already exists")
    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    print(f"Output Zarr directory will be created at '{args.output}'.")

    #inpath = os.path.join(args.xoa_dir, "outs")
    inpath = args.xoa_dir
    if not os.path.exists(inpath):
        raise FileNotFoundError(f"Xenium Onboard Analysis directory '{inpath}' not found")

    spdio.xenium(
        path=inpath,
        cells_boundaries=True,
        nucleus_boundaries=True,
        cells_as_circles=False,
        cells_labels=True,  # recommended
        nucleus_labels=True,  # recommended
        transcripts=True,
        morphology_mip=True,
        morphology_focus=True,
        aligned_images=True,
        cells_table=args.import_adata,  # optional, can be generated from cell boundaries
        gex_only=False
    ).write(args.output)

    print(f"{GREEN}Zarr dataset created successfully.{RESET}")


if __name__ == "__main__":
    main()
