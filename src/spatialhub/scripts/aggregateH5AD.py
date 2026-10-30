import os
import spatialdata as spd

from argparse import ArgumentParser as AP
from pathlib import Path
from spatialhub.scripts.utils import die, print_arguments, RESET, GREEN


def main():
    p = AP(description="Aggregate SpatialData Zarr stack into a raw-counts AnnData using the specified Shapes.")
    p.add_argument("h5adout", help="Output h5ad file path.")
    p.add_argument("--in-zarr", required=True, help="Input Zarr file path.")
    p.add_argument("--points-from", required=True, help="Points to be aggregated.")
    p.add_argument("--shapes-by", required=True, help="Shapes to aggregate by.")
    p.add_argument("--agg-func", default="sum", choices=["mean", "sum", "median"], help="Aggregation function.")
    p.add_argument("--coords", default="global", help="Target coordinate system for the aggregated AnnData.")
    args = p.parse_args()
    print_arguments(args)

    output = Path(args.h5adout)
    if output.exists():
        raise FileExistsError(output)
    
    sdata = spd.read_zarr(args.in_zarr)
    adata = sdata.aggregate(
        values=args.points_from,
        by=args.shapes_by,
        agg_func=args.agg_func,
        target_coordinate_system=args.coords,
        deepcopy=True,
    )["table"]
    
    # Add minimal essential metadata to adata.obs
    # Centroids must be expressed in the same coordinate system as aggregation.
    shapes = sdata.transform_element_to_coordinate_system(
        element_name=args.shapes_by,
        target_coordinate_system=args.coords,
    )
    centroids = shapes.loc[adata.obs_names, "geometry"].centroid
    adata.obs[["array_col", "array_row"]] = [[c.x, c.y] for c in centroids]
    adata.obs["sample_id"] = Path(args.in_zarr).stem
    
    # Store provenance of the aggregated counts matrix in adata.uns
    # to avoid loss of information in case of file name changes.
    adata.uns["counts_mtx_source"] = {
        "input_zarr": str(args.in_zarr),
        "sample_id": Path(args.in_zarr).stem,
        "region_type": "segment",
        "shapes_by": args.shapes_by,
        "points_from": [args.points_from],
        "coords_system": args.coords,
    }

    # Write final AnnData H5AD file
    output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output, compression="gzip")
    
    print(f"{GREEN}Successfully aggregated AnnData.{RESET}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        die("Failed to aggregate AnnData.")
