import geopandas as gpd
import spatialdata as spd

from argparse import ArgumentParser as AP
from spatialhub.scripts.utils import die, print_arguments, RESET, GREEN


def main():
    p = AP(description="Add GeoJSON shapes to the Zarr object.")
    p.add_argument("zarr_path", help="Path to the Zarr directory.")
    p.add_argument("--geojson-path", required=True, help="Path to the GeoJSON file.")
    p.add_argument("--key", required=True, help="Zarr key for the shapes data.")
    p.add_argument("--index", required=True, help="Index column in the GeoJSON table.")
    args = p.parse_args()
    print_arguments(args)

    sdata = spd.read_zarr(args.zarr_path)
    if args.key in sdata.shapes:
        raise KeyError(f"Shape key '{args.key}' already exists.")

    gdf = gpd.read_file(args.geojson_path)
    if gdf.empty:
        raise ValueError("GeoJSON contains no shapes.")
    if gdf.geometry.isna().any() or gdf.geometry.is_empty.any():
        raise ValueError("GeoJSON contains missing or empty geometries.")
    if not gdf.geometry.is_valid.all():
        raise ValueError("GeoJSON contains invalid geometries.")

    gdf.set_crs(None, allow_override=True, inplace=True)

    if args.index not in gdf.columns:
        raise KeyError(f"Index column '{args.index}' not found in GeoJSON.")
    if gdf[args.index].isna().any():
        raise ValueError(f"Missing values in index column '{args.index}'.")
    gdf.set_index(args.index, inplace=True)
    gdf.index = gdf.index.astype(str)
    if not gdf.index.is_unique:
        raise ValueError(f"Duplicate values in index column '{args.index}'.")

    sdata[args.key] = spd.models.ShapesModel.parse(gdf)
    sdata.write_element(args.key)

    print(f"{GREEN}Successfully added shape {args.key}.{RESET}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        die("Failed to add shape.")
