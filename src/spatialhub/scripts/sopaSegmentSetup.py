
#!/usr/bin/env python3

"""
Create a reduced SpatialData Zarr store compatible with SOPA segmentation.

Supports image-based segmentation and transcript-based segmentation,
optionally including an existing segmentation as a prior.
"""

import argparse
from pathlib import Path

import spatialdata as sd
from spatialhub.scripts.utils import die, GREEN, RESET, print_arguments


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input-zarr",
        required=True,
        type=Path,
        help="Path to the input SpatialData Zarr store.",
    )
    parser.add_argument(
        "--output-zarr",
        required=True,
        type=Path,
        help="Path to the reduced SpatialData Zarr store.",
    )
    parser.add_argument(
        "--segmentation-mode",
        choices=["image", "transcripts"],
        default="image",
        help="Type of segmentation to prepare for.",
    )
    parser.add_argument(
        "--image",
        default="image",
        help="Name of the input image element.",
    )
    parser.add_argument(
        "--points",
        default="main",
        help="Name of the input transcripts element.",
    )
    parser.add_argument(
        "--prior",
        default=None,
        help="Name of an existing segmentation to retain.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite an existing output Zarr store.",
    )

    args = parser.parse_args()
    print_arguments(args)


    # ----------------------- Setup ----------------------- #

    inzarr = args.input_zarr
    outzarr = args.output_zarr

    if not inzarr.is_dir():
        raise FileNotFoundError(f"Input Zarr not found: {inzarr}")

    if inzarr.resolve() == outzarr.resolve():
        raise ValueError("Input and output Zarr paths must differ.")

    if outzarr.exists() and not args.overwrite:
        raise FileExistsError(
            f"Output already exists: {outzarr}. "
            "Use --overwrite to replace it."
        )


    # ----------------------- Load SpatialData ----------------------- #

    print(f"Reading SpatialData: {inzarr}")
    sdata = sd.read_zarr(inzarr)
    print(sdata)


    # ----------------------- Select elements ----------------------- #

    elements = []

    if args.image not in sdata.images:
        raise KeyError(
            f"Image element '{args.image}' not found."
        )
    elements.append(args.image)


    if args.segmentation_mode == "transcripts":

        if args.points not in sdata.points:
            raise KeyError(
                f"Transcripts element '{args.points}' not found."
            )
        elements.append(args.points)

        if args.prior is not None:
            if args.prior not in sdata.shapes:
                raise KeyError(
                    f"Prior segmentation '{args.prior}' not found."
                )
        elements.append(args.prior)


    # ----------------------- Subset SpatialData ----------------------- #

    print(f"Retaining elements: {elements}")
    sdata = sdata.subset(elements)

    if args.segmentation_mode == "transcripts":
        if args.points != "transcripts":
            sdata["transcripts"] = sdata[args.points]
            del sdata[args.points]

    print(sdata)


    # ----------------------- Write output ----------------------- #

    outzarr.parent.mkdir(parents=True, exist_ok=True)
    print(f"Writing SOPA-compatible SpatialData: {outzarr}")
    sdata.write(outzarr, overwrite=args.overwrite)

    print("SOPA setup completed.")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        die("Failed to prepare Zarr stack for SOPA segmentation.")