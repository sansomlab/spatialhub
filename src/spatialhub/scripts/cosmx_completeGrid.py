"""
Complete the CosMx field grid with symbolic links to FOV tiles and mock FOVs.

This utility reads a CSV file containing CosMx FOV positions and creates
a symbolic-link-based directory structure for downstream analysis,
filling missing FOVs with a blank mock FOV tile to form a complete grid.

Outputs
-------
- A directory containing symbolic links to the FOV tiles.
- A CSV file describing the complete grid of fields, with the following columns:
    - grid_index: Index of the field in the complete grid (0-based).
    - x_global_px: X-coordinate of the field in global pixel coordinates.
    - y_global_px: Y-coordinate of the field in global pixel coordinates.
    - col_index: Column index of the field in the grid (0-based).
    - row_index: Row index of the field in the grid (0-based).
    - FOV: Original CosMx FOV identifier, or `-1` for mock fields.
"""

import os
import pandas as pd

from argparse import ArgumentParser as AP
from tifffile import imread
from itertools import product
from spatialhub.scripts.utils import print_arguments, RESET, GREEN, die, parse_fov_list


def main():
    p = AP(description="Organise CosMx FOVs into a complete grid.")
    p.add_argument("outdir", help="Output directory for the field symbolic links.")
    p.add_argument("--fov-csv", required=True, help="CosMx FOV position file.")
    p.add_argument("--fov-lst", required=True, help="Comma-separated FOVs to include.")
    p.add_argument("--m2d-pfx", required=True, help="Prefix of Morphology2D paths.")
    p.add_argument("--mock-tiff", required=True, help="Path to a blank mock FOV tile.")
    args = p.parse_args()
    print_arguments(args)

    out_csv = os.path.join(args.outdir, "grid_positions.csv")
    if os.path.exists(out_csv):
        raise FileExistsError(f"file '{out_csv}' already exists")
    os.makedirs(os.path.join(args.outdir, "field_links"), exist_ok=False)

    # Read mock FOV tile to get width and height
    _, height, width = imread(args.mock_tiff).shape
    print(f"Parsed grid dimensions: {width} x {height} pixels.")

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
    fovpos = fovpos.loc[fovs, ["x_global_px", "y_global_px"]].copy()

    # Compute FOV grid coordinates for each FOV, and assign to fovpos
    xmin_px, ymax_px = fovpos["x_global_px"].min(), fovpos["y_global_px"].max()
    ## x, col, width
    fovpos["icol"] = ((fovpos["x_global_px"] - xmin_px) / width + 0.5).astype(int)
    ## y, row, height
    fovpos["irow"] = ((ymax_px - fovpos["y_global_px"]) / height + 0.5).astype(int)

    # Compute grid dimensions and update the field positions to include missing fields
    ncols, nrows = int(fovpos["icol"].max() + 1), int(fovpos["irow"].max() + 1)
    grid2fov = fovpos.reset_index().set_index(["icol", "irow"])["FOV"].to_dict()
    rows = []
    for igrid, (irow, icol) in enumerate(product(range(nrows), range(ncols))):
        if (icol, irow) in grid2fov:
            fov = grid2fov[(icol, irow)]
            x_px, y_px = fovpos.loc[fov, ["x_global_px", "y_global_px"]].tolist()
            src = f"{args.m2d_pfx}{fov:05}.TIF"
        else:
            fov = -1
            src = args.mock_tiff
            x_px, y_px = xmin_px + icol * width, ymax_px - irow * height
        dest = os.path.join(args.outdir, "field_links", f"F{igrid:05}.TIF")
        rows.append(
            {
                "grid_index": igrid,
                "x_global_px": x_px,
                "y_global_px": y_px,
                "col_index": icol,
                "row_index": irow,
                "FOV": fov,
            }
        )
        if not os.path.isfile(src):
            raise FileNotFoundError(f"source file '{src}' does not exist")
        os.symlink(src, dest)
    field_pos = pd.DataFrame(rows)
    field_pos.to_csv(out_csv, index=False)

    print(f"{GREEN}Created {field_pos.shape[0]} field links in {args.outdir}.{RESET}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        die("Failed to organise fields.")
