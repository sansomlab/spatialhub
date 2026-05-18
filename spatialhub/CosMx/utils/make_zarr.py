import argparse
import os
import numpy as np
import pandas as pd
import tifffile
import geopandas as gpd
import spatialdata as spd

from shapely.geometry import Polygon

FOV_WIDTH, FOV_HEIGHT = 4256, 4256

if __name__ == "__main__":
    print("Parsing arguments")
    parser = argparse.ArgumentParser(
        description="Stitch image tiles into a single mosaic."
    )
    parser.add_argument(
        "--out-zarr", required=True, type=str, help="Path to the output Zarr file."
    )
    parser.add_argument(
        "--fovpos-csv",
        required=True,
        help="Path to the CSV file containing FOV positions.",
    )
    parser.add_argument(
        "--tx-csv",
        required=True,
        help="Path to the CSV file containing transcript information.",
    )
    parser.add_argument(
        "--pg-csv",
        default=None,
        help="Path to the CSV file containing polygon annotations.",
    )
    parser.add_argument(
        "--ctrl-probe-regex",
        required=True,
        help="Regex pattern to identify control probes in the transcript table.",
    )
    parser.add_argument(
        "--morph-pfx", type=str, help="Prefix for the morphological image files."
    )
    parser.add_argument(
        "--fov-lst",
        required=True,
        help="Comma-separated FOVs for the object (e.g. '1,2,3').",
    )
    parser.add_argument(
        "--channel-lst",
        required=True,
        help="Comma-separated list of channels to include.",
    )
    parser.add_argument(
        "--fov-width", default=FOV_WIDTH, help="Width of each FOV in pixels."
    )
    parser.add_argument(
        "--fov-height", default=FOV_HEIGHT, help="Height of each FOV in pixels."
    )
    parser.add_argument(
        "--scale-factors", default=None, help="Comma-separated list of scale factors."
    )
    args = parser.parse_args()

    if os.path.exists(args.out_zarr):
        raise FileExistsError(
            f"Output file {args.out_zarr} already exists. Exiting to avoid overwriting."
        )

    channels = [ch.strip() for ch in args.channel_lst.split(",")]
    print(f"Channels to include: {channels}")
    n_channels = len(channels)

    fov_lst = [int(fov.strip()) for fov in args.fov_lst.split(",")]

    if args.scale_factors is not None:
        scale_factors = [int(sf.strip()) for sf in args.scale_factors.split(",")]
    else:
        scale_factors = None

    fov_pos = pd.read_csv(args.fovpos_csv, index_col=0, header=0).loc[fov_lst].copy()

    xmin, ymin = fov_pos[["x_global_px", "y_global_px"]].min(axis=0)
    ymin -= int(args.fov_height)  # y = 0 at the top edge of the bottom row of FOVs
    xmax, ymax = fov_pos[["x_global_px", "y_global_px"]].max(axis=0)
    xmax += int(args.fov_width)  # x = xmax at the left edge of rightmost column of FOVs
    W = xmax - xmin
    H = ymax - ymin
    print(f"Calculated mosaic dimensions: {W} x {H}")

    tile_dict = {}
    for fov in fov_lst:
        path2tif = f"{args.morph_pfx}{fov:05}.TIF"
        tile_dict[fov] = tifffile.imread(path2tif)

    canvas = np.zeros((n_channels, H, W), dtype=tile_dict[fov_lst[0]].dtype)
    for fov, (x, y) in fov_pos[["x_global_px", "y_global_px"]].iterrows():
        tile = tile_dict[fov]
        x_start = x - xmin  # x = xstart at the left edge of the FOV
        y_end = y - ymin  # y = yend at the top edge of the FOV
        x_end = x_start + int(args.fov_width)
        y_start = y_end - int(args.fov_height)
        canvas[:, y_start:y_end, x_start:x_end] = tile[:, ::-1, :]
    images = {
        "image": spd.models.Image2DModel.parse(
            canvas,
            c_coords=channels,
            scale_factors=scale_factors,
        )
    }

    tx_all = pd.read_csv(args.tx_csv)
    tx_all = tx_all.loc[tx_all["fov"].isin(fov_lst)].copy()
    tx_all["x_global_px"] -= xmin
    tx_all["y_global_px"] -= ymin
    idx_ctrl = tx_all["target"].str.contains(args.ctrl_probe_regex, regex=True)
    points = {
        "tx_ctrl": spd.models.PointsModel.parse(
            tx_all[idx_ctrl],
            feature_key="target",
            coordinates={"x": "x_global_px", "y": "y_global_px"},
        ),
        "tx_main": spd.models.PointsModel.parse(
            tx_all[~idx_ctrl],
            feature_key="target",
            coordinates={"x": "x_global_px", "y": "y_global_px"},
        ),
    }

    if args.pg_csv is not None:
        df = pd.read_csv(args.pg_csv)
        df = df.loc[df["fov"].isin(fov_lst)].copy()
        df["x_global_px"] -= xmin
        df["y_global_px"] -= ymin
        segs = []
        for cell, sub in df.groupby("cell"):
            coords = list(zip(sub["x_global_px"], sub["y_global_px"]))
            if len(coords) > 2:
                poly = Polygon(coords)
                segs.append({"cell_id": cell, "geometry": poly})
        shapes = {
            "shapes": spd.models.ShapesModel.parse(gpd.GeoDataFrame(segs, crs=None))
        }
    else:
        shapes = None

    spd.SpatialData(images=images, points=points, shapes=shapes).write(args.out_zarr)
