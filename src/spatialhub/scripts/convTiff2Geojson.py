import numpy as np
import rasterio
import geopandas as gpd

from rasterio.features import shapes
from shapely.geometry import shape
from argparse import ArgumentParser as AP
from spatialhub.scripts.utils import die, print_arguments, RESET, GREEN


def main():
    p = AP(description="Convert a segmentation mask TIFF file to GeoJSON format.")
    p.add_argument("geojson_out", help="Path to the output GeoJSON file.")
    p.add_argument("--in-tiff", required=True, help="Path to the input TIFF file.")
    args = p.parse_args()
    print_arguments(args)

    with rasterio.open(args.in_tiff) as src:
        if src.count != 1:
            raise ValueError("Segmentation TIFF must have a single band.")
        msk = src.read(1)
    if not np.issubdtype(msk.dtype, np.integer):
        raise ValueError(f"Integer mask needed, got {msk.dtype}")
    if msk.max() > np.iinfo(np.int32).max:
        raise ValueError("Segment labels exceed int32 range.")
    if msk.max() == 0:
        raise ValueError("Segment labels must contain positive values.")
    if msk.min() < 0:
        raise ValueError("Segment labels must be non-negative.")
    if not np.any(msk == 0):
        raise ValueError("Segment labels must contain 0 as background label.")

    segments = []
    pg_lst = shapes(msk.astype(np.int32, copy=False), mask=(msk != 0))
    for geom, value in pg_lst:
        segments.append({"label": int(value), "geometry": shape(geom)})
    gdf = gpd.GeoDataFrame(segments, crs=None)
    gdf = gdf.dissolve(by="label")
    gdf["area"] = gdf.area

    _, npx = np.unique(msk, return_counts=True)
    if not np.allclose(gdf.sort_index()["area"].values, npx[1:], atol=1e-9):
        raise ValueError("Polygon areas do not match pixel counts.")

    gdf.to_file(args.geojson_out, driver="GeoJSON")

    print(f"{GREEN}Successfully converted segmentation mask to GeoJSON.{RESET}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        die("Failed to convert segmentation mask to GeoJSON.")
