"""
Stitch selected CosMx FOVs with ASHLAR using acquisition coordinates.

Missing FOVs and irregular tile layouts are supported without
reconstructing a rectangular grid or creating mock tiles.

Outputs
-------
- image.ome.tiff: stitched image.
- FOV.positions.csv: stitched positions, shifts and original coordinates.

This script is adapted from https://gist.github.com/jmuhlich/59e114243278cf27e7f0720c2c956484.
"""

from argparse import ArgumentParser as AP
import os
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from ashlar.fileseries import FileSeriesReader
from ashlar.reg import EdgeAligner, Mosaic, PyramidWriter
from spatialhub.scripts.utils import (
    GREEN, RESET, print_arguments, die, parse_fov_list
)


def main():
    p = AP(description="Stitch specified CosMx FOV tiles using ASHLAR.")
    p.add_argument("outdir", help="Directory for the output files.")
    p.add_argument("--fov-csv", required=True, help="CosMx FOV position file.")
    p.add_argument("--fov-lst", required=True, help="Comma-separated FOVs to include.")
    p.add_argument("--m2d-pfx", required=True, help="Prefix of Morphology2D paths.")
    p.add_argument("--px-size", type=float, default=0.120280945, help="Pixel size.")
    p.add_argument("--overlap", type=float, default=0.01, help="FOV overlap fraction.")
    p.add_argument("--max-shift", type=int, default=15, help="Maximum shift (µm).")
    p.add_argument("--ch-align", type=int, default=0, help="Channel for alignment.")
    p.add_argument("--ch-out", default=None, help="Comma-separated channels to output.")
    args = p.parse_args()
    print_arguments(args)

    # Setup output file and validate arguments
    out_csv = os.path.join(args.outdir, "FOV.positions.csv")
    out_tiff = os.path.join(args.outdir, "image.ome.tiff")
    if os.path.exists(out_csv) or os.path.exists(out_tiff):
        raise FileExistsError(f"output file '{out_csv}' or '{out_tiff}' already exists")
    print(f"Output files will be written to '{out_csv}' and '{out_tiff}'.")
    os.makedirs(args.outdir, exist_ok=True)

    if not 0 <= args.overlap < 1:
        raise ValueError("--overlap must be between 0 and 1")
    if not np.isfinite(args.px_size) or args.px_size <= 0:
        raise ValueError("--px-size must be positive and finite")
    if args.max_shift < 0:
        raise ValueError("--max-shift must be non-negative")


    # Read CosMx FOV positions and filter to requested FOVs
    fovs = parse_fov_list(args.fov_lst)  # Do we want to sort them?! sorted(set(parse_fov_list(args.fov_lst)))
    fovpos = pd.read_csv(args.fov_csv)
    if "fov" in fovpos.columns:
        fovpos.rename(columns={"fov": "FOV"}, inplace=True)

    required = {"FOV", "x_global_px", "y_global_px"}
    missing = required - set(fovpos.columns)
    if missing:
        raise ValueError(f"missing variables in FOV table: {sorted(missing)}")

    fovpos.set_index("FOV", inplace=True)
    if not fovpos.index.is_unique:
        raise ValueError("duplicate FOV identifiers")
    
    missing = sorted(set(fovs) - set(fovpos.index))
    if missing:
        raise ValueError(f"missing FOVs in FOV table: {', '.join(map(str, missing))}")
    fovpos = fovpos.loc[fovs, ["x_global_px", "y_global_px"]].copy()

    if not np.isfinite(fovpos.to_numpy(dtype=float)).all():
        raise ValueError("invalid CosMx FOV coordinates")


    # Parse channels to write
    ch_lst = list(map(int, args.ch_out.split(","))) if args.ch_out is not None else None


    # Create temporary symlinks to tiles (no mock FOV needed)
    with TemporaryDirectory(prefix="ashlar_fovs_") as field_dir:
        for fov in fovs:

            src = os.path.abspath(f"{args.m2d_pfx}{fov:05}.TIF")
            if not os.path.isfile(src):
                raise FileNotFoundError(f"Source file '{src}' does not exist")

            dest = os.path.join(field_dir, f"F{fov:05}.TIF")
            os.symlink(src, dest)

        # Read the field tiles using FileSeriesReader (will be overridden)
        reader = FileSeriesReader(
            field_dir,
            pattern="F{series}.TIF",
            width=len(fovs),
            height=1,
            overlap=args.overlap,
            pixel_size=args.px_size,
        )

        # Save original FOV dimensions (height, width) in pixels - useful for cosmx_makeZarr
        height, width = map(int, reader.metadata.size)
        out_size = os.path.join(args.outdir, "FOV.size.px.txt")
        with open(out_size, "w") as f:
            f.write(f"{width},{height}\n")
        print(f"FOV dimensions: {width} x {height} pixels")

        # Match coordinates to ASHLAR's internal image ordering
        series_fovs = [int(s) for _, s in reader.metadata.all_series]
        if sorted(series_fovs) != sorted(fovs):
            raise ValueError("ASHLAR FileSeriesReader FOV identifiers do not match input FOVs")

        pos = fovpos.loc[series_fovs]

        # Convert CosMx coordinates into a common coordinate system
        # with X increasing rightwards and Y increasing downwards
        # and applying the fixed overlap correction
        # (WARNING: may need changing if CosMx changes its FOV coordinate system)
        x = pos["x_global_px"].to_numpy(dtype=float)
        y = pos["y_global_px"].to_numpy(dtype=float)
        scale = 1.0 - args.overlap
        positions = np.column_stack([((y.max() - y)) * scale, (x - x.min()) * scale])

        # Replace FileSeriesReader's regular-grid positions with actual CosMx positions
        if positions.shape != reader.metadata.positions.shape:
            raise ValueError("Position array shape mismatch")
        reader.metadata._positions = positions
        if not np.allclose(reader.metadata.positions, positions):
            raise ValueError("Failed to update ASHLAR positions")
        

        # Align the field tiles using EdgeAligner
        aligner = EdgeAligner(reader, channel=args.ch_align, max_shift=args.max_shift)
        
        print(f"Stitching {len(fovs)} FOVs")
        print(
            "Candidate neighbour pairs:",
            aligner.neighbors_graph.number_of_edges()
        )

        # Run the alignment process
        aligner.run()

        accepted = sum(
            np.isfinite(error)
            for _, error in aligner._cache.values()
        )
        print(f"Accepted registrations: {accepted}")

        # Write out the stitched mosaic to an OME-TIFF file
        mosaic = Mosaic(aligner, aligner.mosaic_shape, channels=ch_lst)
        writer = PyramidWriter([mosaic], out_tiff)
        writer.run()
        print(f"Wrote output to '{out_tiff}'.")

        # Write out the positions and shifts to a CSV file
        coords_df = pd.DataFrame(
            np.hstack([aligner.positions, aligner.shifts]),
            columns=["Position_Y", "Position_X", "Shift_Y", "Shift_X"],
            index=pd.Index(series_fovs, name="FOV"),
        ).round(5)

        coords_df["x_global_px"] = pos["x_global_px"]
        coords_df["y_global_px"] = pos["y_global_px"]

        coords_df = coords_df[
            [
                "Position_X", "Position_Y",
                "Shift_X", "Shift_Y",
                "x_global_px", "y_global_px",
            ]
        ]

        coords_df.to_csv(out_csv)
        print(f"Wrote output to '{out_csv}'.")

    print(f"{GREEN}ASHLAR stitching completed successfully.{RESET}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        die("Failed to run ASHLAR.")
