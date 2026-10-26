import os
import spatialdata as spd

from argparse import ArgumentParser as AP
from pathlib import Path
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
    
    # Add key metadata to adata.obs
    centroids = sdata[args.shapes_by].loc[adata.obs_names, "geometry"].centroid
    adata.obs[["array_col", "array_row"]] = [[c.x, c.y] for c in centroids]
    
    sample_id = Path(args.in_zarr).stem
    adata.obs["sample_id"] = sample_id
    
    adata.layers["counts"] = adata.X.copy()  # Ensure counts layer is present
    
    # Store provenance of the aggrgated counts matrix in adata.uns
    # to avoid loss of information in case of file name changes.
    adata.uns["counts_mtx_source"] = {
        "input_zarr": str(args.in_zarr),
        "sample_id": sample_id,
        "points_from": [
            point.strip()
            for point in args.points_from.split(",")
            if point.strip()
        ],
        "shapes_by": args.shapes_by,
        "agg_func": args.agg_func,
        "coords_system": args.coords,
    }

    # Write final AnnData H5AD file
    adata.write(args.h5adout, compression="gzip")

    print(f"{GREEN}Successfully aggregated AnnData.{RESET}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        die("Failed to aggregate AnnData.")
