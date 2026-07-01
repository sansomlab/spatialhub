import os
import pandas as pd

from argparse import ArgumentParser as AP
from tifffile import imread, imwrite
from itertools import product
from spatialhub.scripts.utils import print_arguments, RED, RESET, GREEN, die


def main():
    p = AP(description="Create a symbolic-link-based directory for a CosMx sample.")
    # Required arguments
    p.add_argument("outdir", help="Output directory for the reorganised structure.")
    p.add_argument("--fov-csv", required=True, help="CosMx FOV position file.")
    p.add_argument("--fov-lst", required=True, help="Comma-separated FOVs to include.")
    p.add_argument("--m2d-pfx", required=True, help="Prefix of Morphology2D paths.")
    p.add_argument("--mock-path", required=True, help="Path to a blank mock FOV tile.")
    args = p.parse_args()
    print_arguments(args)

    if os.path.exists(args.outdir):
        raise FileExistsError(f"directory '{args.outdir}' already exists")
    os.makedirs(os.path.join(args.outdir, "FOVs"), exist_ok=False)

    # Read mock FOV tile to get width and height
    height, width = imread(args.mock_path).shape
    print(f"Parsed FOV tile dimensions: {width} x {height} pixels.")

    # Parse FOVs to include
    fovs = []
    for fov_str in args.fov_lst.split(","):
        if "-" in fov_str:
            try:
                start, end = map(int, fov_str.split("-"))
            except ValueError:
                raise ValueError(f"invalid FOV range '{fov_str}'")
            if start > end:
                raise ValueError(f"invalid FOV range '{fov_str}': start > end")
            fovs.extend(range(start, end + 1))
        else:
            try:
                fovs.append(int(fov_str))
            except ValueError:
                raise ValueError(f"invalid FOV '{fov_str}'")

    # Read FOV positions and filter to requested FOVs
    fovpos = pd.read_csv(args.fov_csv)
    if not fovpos.columns.str.contains(r"FOV|fov").any():
        raise ValueError("missing 'FOV' or 'fov' columns in FOV table")
    fovpos.rename(columns={"FOV": "fov"}, inplace=True)
    fovpos.set_index("fov", inplace=True)

    if not "x_global_px" in fovpos.columns or not "y_global_px" in fovpos.columns:
        raise ValueError("missing 'x_global_px' and 'y_global_px' columns in FOV table")
    fovpos = fovpos[["x_global_px", "y_global_px"]]
    fovpos = fovpos.loc[fovs].copy()

    # Compute FOV grid coordinates for each FOV, and assign to fovpos
    xmin_px, ymax_px = fovpos["x_global_px"].min(), fovpos["y_global_px"].max()
    ## x, col, width
    fovpos["icol"] = ((fovpos["x_global_px"] - xmin_px) / width + 0.5).astype(int)
    ## y, row, height
    fovpos["irow"] = ((ymax_px - fovpos["y_global_px"]) / height + 0.5).astype(int)

    # Compute grid dimensions and make a table for symbolic link creation
    ncols, nrows = fovpos["icol"].max() + 1, fovpos["irow"].max() + 1
    grid2fov = fovpos.reset_index().set_index(["icol", "irow"])["fov"].to_dict()
    rows = []
    for ifov, (irow, icol) in enumerate(product(range(nrows), range(ncols))):
        if (icol, irow) in grid2fov:
            src = f"{args.m2d_pfx}{grid2fov[(icol, irow)]:05}.TIF"
        else:
            src = args.mock_path
        dest = os.path.join(args.outdir, "FOVs", f"F{ifov:05}.TIF")
        rows.append({"icol": icol, "irow": irow, "src": src, "dest": dest})
        os.symlink(src, dest)
    src2dest = pd.DataFrame(rows)
    src2dest.to_csv(os.path.join(args.outdir, "fov_links.summary.csv"), index=False)
    print(f"{GREEN}Created links for {len(src2dest)} FOVs in {args.outdir}.{RESET}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        die("Failed to organise FOVs.")
