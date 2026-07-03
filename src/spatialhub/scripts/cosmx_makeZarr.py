import os
import pandas as pd
import geopandas as gpd
import tifffile
import spatialdata as spd

from argparse import ArgumentParser as AP
from shapely import Polygon
from spatialhub.scripts.utils import die, print_arguments, RESET, YELLOW


def transform_coordinates(points_df, fov2pos):
    points_df["x"] = points_df.apply(
        lambda r: fov2pos[r["FOV"]]["Position_X"]
        + r["x_global_px"]
        - fov2pos[r["FOV"]]["x_global_px"],
        axis=1,
    )
    points_df["y"] = points_df.apply(
        lambda r: fov2pos[r["FOV"]]["Position_Y"]
        + fov2pos[r["FOV"]]["y_global_px"]
        - r["y_global_px"],
        axis=1,
    )
    return points_df


def main():
    p = AP(description="Prepare a Zarr object from CosMx data.")
    # Required arguments
    p.add_argument("--out-zarr", required=True, help="Output Zarr directory.")
    p.add_argument("--fov-csv", required=True, help="CosMx FOV position CSV file.")
    p.add_argument("--tx-csv", required=True, help="CosMx transcript CSV file.")
    p.add_argument("--morph-tiff", required=True, help="Assembled/stitched TIFF image.")
    # Optional arguments
    p.add_argument("--pg-csv", default=None, help="CosMx polygon CSV file.")
    p.add_argument("--channels", default=None, help="Comma-separated channel names.")
    p.add_argument("--ctrl-regex", default=None, help="Regex for control probes.")
    p.add_argument("--other-regex", default=None, help="Regex for other probes.")
    p.add_argument("--scales", default=None, help="Comma-separated image scales.")
    args = p.parse_args()
    print_arguments(args)

    if os.path.exists(args.out_zarr):
        raise FileExistsError(f"output file {args.out_zarr} already exists")
    out_dir = os.path.dirname(args.out_zarr)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    print(f"Output Zarr directory will be created at '{args.out_zarr}'.")

    if not os.path.exists(args.fov_csv):
        raise FileNotFoundError(f"FOV position file '{args.fov_csv}' not found")
    if not os.path.exists(args.tx_csv):
        raise FileNotFoundError(f"Transcript file '{args.tx_csv}' not found")
    if not os.path.exists(args.morph_tiff):
        raise FileNotFoundError(f"TIFF file '{args.morph_tiff}' not found")

    img = tifffile.imread(args.morph_tiff)
    if args.channels:
        channels = [ch.strip() for ch in args.channels.split(",")]
    else:
        channels = list(range(img.shape[0]))

    sf = list(map(int, args.scales.split(","))) if args.scales else None
    print(f"Image shape: {img.shape}\nChannels: {channels}\nScale factors: {sf}")

    # Create a Image2DModel object for the TIFF image, with the specified channels and scale factors
    images = {
        "image": spd.models.Image2DModel.parse(img, c_coords=channels, scale_factors=sf)
    }

    # Read the FOV position CSV
    fovpos = pd.read_csv(args.fov_csv, index_col="FOV", header=0)
    missing = {"Position_X", "Position_Y"} - set(fovpos.columns)
    if missing:
        raise ValueError(f"missing {', '.join(sorted(missing))} in FOV table")
    fov2pos = fovpos[
        ["Position_X", "Position_Y", "x_global_px", "y_global_px"]
    ].to_dict(orient="index")

    # Read the transcript CSV and compute global coordinates for each transcript
    tx = pd.read_csv(args.tx_csv, header=0)
    if "fov" in tx.columns:
        tx.rename(columns={"fov": "FOV"}, inplace=True)
    missing = {"FOV", "x_global_px", "y_global_px", "target"} - set(tx.columns)
    if missing:
        raise ValueError(f"missing {', '.join(sorted(missing))} in tx table")
    tx = tx.loc[tx["FOV"].isin(fov2pos.keys())].copy()
    tx = transform_coordinates(tx, fov2pos)

    # Split the transcripts into main, control, and other categories based on regex patterns
    idx_ctrl = (
        tx["target"].str.contains(args.ctrl_regex, na=False)
        if args.ctrl_regex
        else pd.Series(False, index=tx.index)
    )
    idx_other = (
        tx["target"].str.contains(args.other_regex, na=False)
        if args.other_regex
        else pd.Series(False, index=tx.index)
    )
    idx_main = ~idx_ctrl & ~idx_other

    # Create PointsModel objects for each category of transcripts
    points = {
        "main": spd.models.PointsModel.parse(
            tx.loc[idx_main], feature_key="target", coordinates={"x": "x", "y": "y"}
        ),
        "ctrl": spd.models.PointsModel.parse(
            tx.loc[idx_ctrl], feature_key="target", coordinates={"x": "x", "y": "y"}
        ),
        "other": spd.models.PointsModel.parse(
            tx.loc[idx_other], feature_key="target", coordinates={"x": "x", "y": "y"}
        ),
    }

    # If a polygon CSV is provided, read it and create a ShapesModel object
    if args.pg_csv is not None:
        pg = pd.read_csv(args.pg_csv, header=0)
        if "fov" in pg.columns:
            pg.rename(columns={"fov": "FOV"}, inplace=True)
        missing = {"FOV", "cell", "x_global_px", "y_global_px"} - set(pg.columns)
        if missing:
            raise ValueError(f"missing {', '.join(sorted(missing))} in pg table")
        pg = pg.loc[pg["FOV"].isin(fov2pos.keys())].copy()
        pg = transform_coordinates(pg, fov2pos)
        segments = []
        for cell, sub in pg.groupby("cell"):
            coords = list(zip(sub["x"], sub["y"]))
            if len(coords) > 2:
                poly = Polygon(coords)
                segments.append({"cell_id": cell, "geometry": poly})
        nsk = pg["cell"].nunique() - len(segments)
        if nsk > 0:
            print(f"{YELLOW}Warning: {nsk} segments have insufficient points.{RESET}")
        shapes = {
            "atomx": spd.models.ShapesModel.parse(
                gpd.GeoDataFrame(segments, geometry="geometry", crs=None)
            )
        }
    else:
        shapes = None

    spd.SpatialData(images=images, points=points, shapes=shapes).write(args.out_zarr)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        die("Failed to create Zarr dataset.")
