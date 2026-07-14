import os
import spatialdata as spd

from argparse import ArgumentParser as AP
from spatialhub.scripts.utils import die, print_arguments, RESET, GREEN


def main():
    p = AP(description="Aggregate to AnnData object from SpatialData Zarr file.")
    p.add_argument("h5adout", help="Output h5ad file path.")
    p.add_argument("--in-zarr", help="Input Zarr file path.")
    p.add_argument("--points-from", required=True, help="Points to be aggregated.")
    p.add_argument("--shapes-by", required=True, help="Shapes to aggregate by.")
    p.add_argument(
        "--agg-func",
        default="sum",
        choices=["mean", "sum", "median"],
        help="Aggregation function.",
    )
    p.add_argument(
        "--coords",
        default="global",
        help="Target coordinate system for the aggregated AnnData.",
    )
    args = p.parse_args()
    print_arguments(args)

    if os.path.exists(args.h5adout):
        raise FileExistsError(f"Output file {args.h5adout} already exists.")

    sdata = spd.read_zarr(args.in_zarr)
    adata = sdata.aggregate(
        values=args.points_from,
        by=args.shapes_by,
        agg_func=args.agg_func,
        target_coordinate_system=args.coords,
        deepcopy=True,
    )["table"]
    centroids = sdata[args.shapes_by].loc[adata.obs_names, "geometry"].centroid
    adata.obs[["array_col", "array_row"]] = [[c.x, c.y] for c in centroids]
    adata.layers["counts"] = adata.X.copy()  # Ensure counts layer is present
    adata.write(args.h5adout, compression="gzip")

    print(f"{GREEN}Successfully aggregated AnnData.{RESET}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        die("Failed to aggregate AnnData.")
