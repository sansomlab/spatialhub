import os
import spatialdata_io as spdio

from argparse import ArgumentParser as AP
from spatialhub.scripts.utils import die, print_arguments, RESET, GREEN


def main():
    p = AP(description="Create Zarr dataset from Visium HD data.")
    p.add_argument("output", help="Output Zarr directory.")
    p.add_argument("--sr-dir", required=True, help="Path to the SpaceRanger directory.")
    p.add_argument("--capture-id", required=True, help="Dataset ID (e.g., 'sample1').")
    p.add_argument("--fullres-img", default=None, help="Path to full-resolution image.")
    p.add_argument("--use-raw", action="store_true", help="Use raw counts.")
    args = p.parse_args()
    print_arguments(args)

    if os.path.exists(args.output):
        raise FileExistsError(f"output file {args.output} already exists")
    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    print(f"Output Zarr directory will be created at '{args.output}'.")

    inpath = os.path.join(args.sr_dir, "outs")
    if not os.path.exists(inpath):
        raise FileNotFoundError(f"Space Ranger directory '{inpath}' not found")

    if args.fullres_img and not os.path.exists(args.fullres_img):
        raise FileNotFoundError(f"Full-resolution image '{args.fullres_img}' not found")

    spdio.visium_hd(
        path=inpath,
        dataset_id=args.capture_id,
        filtered_counts_file=not args.use_raw,
        load_segmentations_only=False,
        load_nucleus_segmentations=True,
        bins_as_squares=True,
        annotate_table_by_labels=True,
        fullres_image_file=args.fullres_img,
        load_all_images=True,
        var_names_make_unique=True,
    ).write(args.output)

    print(f"{GREEN}Zarr dataset created successfully.{RESET}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        die("Failed to create Zarr dataset.")
