import os
import spatialdata as spd

from argparse import ArgumentParser as AP
from spatialhub.scripts.utils import die, print_arguments, RESET, GREEN, YELLOW


def main():
    p = AP(description="Aggregate to AnnData object from SpatialData Zarr file.")
    p.add_argument("h5adout", help="Output h5ad file path.")
    p.add_argument("--in-zarr", required=True, help="Input Zarr file path.")
    p.add_argument("--points-from", required=True, help="Points to be aggregated.")
    p.add_argument("--shapes-by", required=True, help="Shapes to aggregate by.")
    p.add_argument("--seg-idx", default=None, help="Segment index column.")
    p.add_argument(
        "--agg-func",
        default="sum",
        choices=["mean", "sum", "median"],
        help="Aggregation function.",
    )
    p.add_argument(
        "--agg-name",
        default="counts",
        help="Layer name for storing the aggregated values.",
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

    if args.seg_idx is not None:
        if sdata[args.shapes_by][args.seg_idx].isna().any():
            raise ValueError(
                f"Segment index column '{args.seg_idx}' contains missing values."
            )
        if not sdata[args.shapes_by][args.seg_idx].is_unique:
            raise ValueError(
                f"Segment index column '{args.seg_idx}' contains duplicate values."
            )
        sdata[args.shapes_by].set_index(args.seg_idx, inplace=True)

    adata = sdata.aggregate(
        values=args.points_from,
        by=args.shapes_by,
        agg_func=args.agg_func,
        target_coordinate_system=args.coords,
        deepcopy=True,
    )["table"]
    shapes_t = spd.transform(sdata[args.shapes_by], to_coordinate_system=args.coords)
    print(
        f"{YELLOW}Segment centroids and area will be computed in coordinate system "
        f"'{args.coords}' (CRS: '{shapes_t.crs}').{RESET}"
    )

    centroids = shapes_t.loc[adata.obs_names, "geometry"].centroid
    adata.obs["centroid.x"] = centroids.x.to_numpy()
    adata.obs["centroid.y"] = centroids.y.to_numpy()
    adata.obs["area"] = shapes_t.loc[adata.obs_names, "geometry"].area.to_numpy()

    if args.agg_name in adata.layers:
        raise ValueError(f"Layer '{args.agg_name}' already exists.")
    adata.layers[args.agg_name] = adata.X.copy()

    adata.write(args.h5adout, compression="gzip")

    print(f"{GREEN}Successfully aggregated AnnData.{RESET}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        die("Failed to aggregate AnnData.")
