#!/usr/bin/env python3

"""
Create a blank CosMx FOV TIFF image.

This utility generates a blank FOV tile filled with zeros to serve as a
placeholder for missing FOVs in a CosMx dataset.

Supported data types:
    - u8
    - u16 (default)
    - u32

Examples:
    create_mock_fov.py blank.tif
    create_mock_fov.py blank.tif --width 2048 --height 2048 --dtype u8
"""

import os
import numpy as np

from argparse import ArgumentParser as AP
from tifffile import imwrite
from spatialhub.scripts.utils import print_arguments, RESET, GREEN, die

S = 4256
DTYPE_MAP = {"u16": np.uint16, "u8": np.uint8, "u32": np.uint32}


def main():
    p = AP(description="Create a mock FOV tile for CosMx.")
    p.add_argument("output", help="Output path for the mock FOV tile.")
    p.add_argument("--width", type=int, default=S, help=f"FOV width in pixels [{S}].")
    p.add_argument("--height", type=int, default=S, help=f"FOV height in pixels [{S}].")
    p.add_argument("--dtype", default="u16", help="FOV data type [u16 | u8 | u32].")
    args = p.parse_args()
    print_arguments(args)

    if os.path.exists(args.output):
        raise FileExistsError(f"file '{args.output}' already exists")

    try:
        fov_dtype = DTYPE_MAP[args.dtype]
    except KeyError:
        raise ValueError(f"unknown dtype '{args.dtype}'")
    print(
        f"Creating blank FOV tile at {args.output}",
        f"with dimensions ({args.width}, {args.height}) and dtype np.{args.dtype.replace('u', 'uint')}.",
    )

    blank = np.zeros((args.height, args.width), dtype=fov_dtype)
    imwrite(args.output, blank)
    print(f"{GREEN}Blank FOV tile created.{RESET}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        die("Failed to create blank FOV tile.")
