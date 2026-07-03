"""
Assembly of FOVs for CosMx data based on a FOV position CSV file.

This utility reads a CSV file containing the positions of CosMx FOVs and generates
a TIFF image with the FOVs arranged according to their specified positions.
It also creates a CSV file describing the transformed coordinates of the FOVs in
the assembled image.

Outputs
-------
- A single assembled image containing all FOVs arranged according to their positions.
- A CSV file describing the transformed coordinates of the FOVs in the assembled image,
  with the following columns:
    - FOV: Original CosMx FOV identifier.
    - x_global_px: X-coordinate of the FOV in global pixel coordinates.
    - y_global_px: Y-coordinate of the FOV in global pixel coordinates.
"""

import os
import numpy as np
import pandas as pd
import tifffile

from argparse import ArgumentParser as AP
from spatialhub.scripts.utils import print_arguments, die, parse_fov_list


def main():
    p = AP(description="Assemble CosMx FOVs based on a position CSV file.")
    p.add_argument("outdir", help="Output directory.")
    p.add_argument("--fov-csv", required=True, help="CosMx FOV position file.")
    p.add_argument("--fov-lst", required=True, help="Comma-separated FOVs to include.")
    p.add_argument("--m2d-pfx", required=True, help="Prefix of Morphology2D paths.")
    p.add_argument("--ch-out", default=None, help="Comma-separated channels to output.")
    args = p.parse_args()
    print_arguments(args)

    out_tiff = os.path.join(args.outdir, "assembled.ome.tiff")
    out_csv = os.path.join(args.outdir, "assembled.positions.csv")
    if os.path.exists(out_tiff) or os.path.exists(out_csv):
        raise FileExistsError(f"Output file '{out_csv}' or '{out_tiff}' already exists")
    os.makedirs(args.outdir, exist_ok=True)
    print(f"Output files will be written to '{out_csv}' and '{out_tiff}'.")

    # Parse FOVs to include
    fovs = parse_fov_list(args.fov_lst)

    # Read FOV positions and filter to requested FOVs
    fovpos = pd.read_csv(args.fov_csv)
    if "fov" in fovpos.columns:
        fovpos.rename(columns={"fov": "FOV"}, inplace=True)
    if "FOV" not in fovpos.columns:
        raise ValueError("missing 'FOV' or 'fov' column in FOV table")
    fovpos.set_index("FOV", inplace=True)

    missing = {"x_global_px", "y_global_px"} - set(fovpos.columns)
    if missing:
        raise ValueError(f"missing columns in FOV table: {', '.join(sorted(missing))}")
    missing = sorted(set(fovs) - set(fovpos.index))
    if missing:
        raise ValueError(f"missing FOVs in FOV table: {', '.join(map(str, missing))}")
    fovpos = fovpos.loc[fovs, ["x_global_px", "y_global_px"]].astype(int).copy()

    # Read the TIF files for each FOV and store them in a dictionary
    edim = None
    tile_dict = {}
    for fov in fovs:
        path2tif = f"{args.m2d_pfx}{fov:05}.TIF"
        tile_dict[fov] = tifffile.imread(path2tif)
        d = tile_dict[fov].shape
        if edim is None:
            edim = d
        elif d != edim:
            raise ValueError(f"FOV{fov} dimension mismatch: got {d}, expected {edim}")
    fov_height, fov_width = int(edim[1]), int(edim[2])

    # Calculate the mosaic dimensions based on the FOV positions
    ## The FOV coordinates are given as the upper-left corner of the FOV,
    ## but the FOVs are arranged from the lower-left corner of the mosaic.
    xmin, ymin = fovpos[["x_global_px", "y_global_px"]].min(axis=0)
    xmax, ymax = fovpos[["x_global_px", "y_global_px"]].max(axis=0)
    W = xmax - xmin + fov_width
    H = ymax - ymin + fov_height
    print(f"Calculated mosaic dimensions: {W} x {H}")

    # Parse channels to write
    ch_lst = (
        list(map(int, args.ch_out.split(",")))
        if args.ch_out is not None
        else range(edim[0])
    )

    rows = []
    canvas = np.zeros((len(ch_lst), H, W), dtype=tile_dict[fovs[0]].dtype)
    for fov, (x, y) in fovpos[["x_global_px", "y_global_px"]].iterrows():
        tile = tile_dict[fov]
        x_start = x - xmin  # x = xstart at the left edge of the FOV
        y_start = ymax - y  # y = ystart at the top edge of the FOV
        x_end = x_start + fov_width
        y_end = y_start + fov_height
        canvas[:, y_start:y_end, x_start:x_end] = tile[ch_lst, :, :]
        rows.append({"FOV": fov, "X_Position": x_start, "Y_Position": y_start})
    tifffile.imwrite(out_tiff, canvas, ome=True, metadata={"axes": "CYX"})
    pd.DataFrame(rows).to_csv(out_csv, index=False)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        die("Failed to assemble FOVs.")
