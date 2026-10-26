import os
import anndata
import spatialdata as spd

from pathlib import Path
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
        if args.seg_idx not in sdata[args.shapes_by].columns:
            raise ValueError(
                f"Segment index column '{args.seg_idx}' not found in shapes."
            )
        if sdata[args.shapes_by][args.seg_idx].isna().any():
            raise ValueError(f"Missing values in segment index '{args.seg_idx}'.")

        sdata[args.shapes_by].set_index(args.seg_idx, inplace=True)
        sdata[args.shapes_by].index = sdata[args.shapes_by].index.astype(str)

        if not sdata[args.shapes_by].index.is_unique:
            raise ValueError(f"Duplicate values in segment index '{args.seg_idx}'.")

    points_from_list = [x.strip() for x in args.points_from.split(",") if x.strip()]
    adata_dict = {}
    for points_from in points_from_list:
        if points_from not in sdata:
            raise KeyError(f"Points '{points_from}' not found in Zarr file.")
        adata_dict[points_from] = sdata.aggregate(
            values=points_from,
            by=args.shapes_by,
            agg_func=args.agg_func,
            target_coordinate_system=args.coords,
            deepcopy=True,
        )["table"]
    adata = anndata.concat(
        adata_dict, axis="var", join="outer", merge="same", label="points_from"
    )
    if not adata.var_names.is_unique:
        raise ValueError("Duplicated feature names after concatenation.")
    if not adata.obs_names.is_unique:
        raise ValueError("Duplicated segment IDs after aggregation.")

    shapes_t = spd.transform(sdata[args.shapes_by], to_coordinate_system=args.coords)
    if not adata.obs_names.isin(shapes_t.index).all():
        raise ValueError("Some aggregated segments are missing from shapes.")
    print(
        f"{YELLOW}Segment centroids and area will be computed in coordinate system "
        f"'{args.coords}' (CRS: '{shapes_t.crs}').{RESET}"
    )
    geometries = shapes_t.loc[adata.obs_names].geometry
    centroids = geometries.centroid
    adata.obs["centroid.x"] = centroids.x.to_numpy()
    adata.obs["centroid.y"] = centroids.y.to_numpy()
    adata.obs["area"] = geometries.area.to_numpy()

    if args.agg_name in adata.layers:
        raise ValueError(f"Layer '{args.agg_name}' already exists.")
    adata.layers[args.agg_name] = adata.X.copy()

    sample_id = Path(args.in_zarr).name.removesuffix(".zarr")
    adata.obs["sample_id"] = sample_id

    adata.uns["counts_mtx_source"] = {
        "input_zarr": str(args.in_zarr),
        "sample_id": sample_id,
        "points_from": points_from_list,
        "shapes_by": args.shapes_by,
        "agg_func": args.agg_func,
        "coords_system": args.coords,
    }

    adata.write(args.h5adout, compression="gzip")

    print(f"{GREEN}Successfully aggregated AnnData.{RESET}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        die("Failed to aggregate AnnData.")
