
#!/usr/bin/env python3

"""
Register a resolved SOPA segmentation into an original SpatialData Zarr.
"""

import argparse
from pathlib import Path

import numpy as np
import spatialdata as sd
from spatialdata.transformations import get_transformation


# ----------------------- Parameters ----------------------- #

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--input-zarr", required=True, type=Path)
    parser.add_argument("--output-zarr", required=True, type=Path)
    parser.add_argument("--source-key", required=True)
    parser.add_argument("--mask-key", required=True)
    parser.add_argument("--image-key", default="image")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace an existing segmentation with the same mask key.",
    )

    args = parser.parse_args()

    source_path = args.input_zarr.resolve()
    destination_path = args.output_zarr.resolve()

    if not source_path.is_dir():
        raise FileNotFoundError(source_path)

    if not destination_path.is_dir():
        raise FileNotFoundError(destination_path)

    if source_path == destination_path:
        raise ValueError("Input and output Zarr must be different.")

    # ----------------------- Load SpatialData ----------------------- #

    print(f"Reading segmentation from: {source_path}")
    sdata_seg = sd.read_zarr(source_path)

    print(f"Reading original SpatialData: {destination_path}")
    sdata_original = sd.read_zarr(destination_path)

    if args.source_key not in sdata_seg.shapes:
        raise KeyError(
            f"Segmentation '{args.source_key}' not found."
        )

    if args.mask_key in sdata_original:
        if not args.overwrite:
            raise FileExistsError(
                f"Element '{args.mask_key}' already exists. "
                "Use --overwrite to replace it."
            )

        if args.mask_key not in sdata_original.shapes:
            raise ValueError(
                f"Element '{args.mask_key}' exists but is not a shapes element."
            )

        print(f"Overwriting existing segmentation: {args.mask_key}")
        sdata_original.delete_element_from_disk(args.mask_key)
        del sdata_original[args.mask_key]

    if 'image' not in sdata_seg.images:
        raise KeyError(
            f"Image 'image' not found in reduced Zarr."
        )

    if args.image_key not in sdata_original.images:
        raise KeyError(
            f"Image '{args.image_key}' not found in original Zarr."
        )

    # ----------------------- Validate spatial compatibility ----------------------- #

    source_image = sdata_seg.images['image']
    original_image = sdata_original.images[args.image_key]

    source_transforms = get_transformation(
        source_image, get_all=True
    )
    original_transforms = get_transformation(
        original_image, get_all=True
    )

    if set(source_transforms) != set(original_transforms):
        raise ValueError(
            "Source and original images have different coordinate systems."
        )

    for cs in source_transforms:
        matrix_source = source_transforms[cs].to_affine_matrix(
            input_axes=("x", "y"),
            output_axes=("x", "y"),
        )

        matrix_original = original_transforms[cs].to_affine_matrix(
            input_axes=("x", "y"),
            output_axes=("x", "y"),
        )

        if not np.allclose(
            matrix_source, matrix_original,
            rtol=1e-6, atol=1e-6
        ):
            raise ValueError(
                f"Image transformations differ in coordinate system '{cs}'."
            )

    # ----------------------- Register mask ----------------------- #

    segmentation = sdata_seg.shapes[args.source_key].copy()

    shape_transforms = get_transformation(
        segmentation, get_all=True
    )

    if not shape_transforms:
        raise ValueError(
            "Segmentation has no spatial transformations."
        )

    if not set(shape_transforms).issubset(original_transforms):
        raise ValueError(
            "Segmentation references coordinate systems "
            "not present in the original image."
        )

    n_cells = len(segmentation)

    print(f"Registering {n_cells} cells as '{args.mask_key}'")

    sdata_original[args.mask_key] = segmentation
    sdata_original.write_element(args.mask_key)

    # ----------------------- Validate output ----------------------- #

    sdata_check = sd.read_zarr(destination_path)

    if args.mask_key not in sdata_check.shapes:
        raise RuntimeError(
            f"Registered segmentation '{args.mask_key}' not found."
        )

    registered = sdata_check.shapes[args.mask_key]

    if len(registered) != n_cells:
        raise RuntimeError(
            "Registered segmentation has an unexpected cell count."
        )

    registered_transforms = get_transformation(
        registered, get_all=True
    )

    if set(shape_transforms) != set(registered_transforms):
        raise RuntimeError(
            "Coordinate systems changed during registration."
        )

    print("Segmentation registration completed successfully.")


if __name__ == "__main__":
    main()
