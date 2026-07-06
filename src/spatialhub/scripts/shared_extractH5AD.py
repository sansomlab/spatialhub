import spatialdata as spd

from argparse import ArgumentParser as AP
from spatialhub.scripts.utils import die, print_arguments, RESET, GREEN


def main():
    p = AP(description="Extract AnnData object from SpatialData Zarr file.")
    p.add_argument("h5adout", help="Output h5ad file path.")
    p.add_argument("--in-zarr", help="Input Zarr file path.")
    p.add_argument("--extract", required=True, help="Table name to extract.")
    args = p.parse_args()
    print_arguments(args)

    sdata = spd.read_zarr(args.in_zarr)

    if args.extract not in sdata:
        raise ValueError(f"Table '{args.extract}' not found in Zarr file.")
    sdata[args.extract].write(args.h5adout, compression="gzip")

    print(f"{GREEN}Successfully extracted AnnData.{RESET}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        die("Failed to extract AnnData.")
