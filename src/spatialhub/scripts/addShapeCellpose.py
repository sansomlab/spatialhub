import spatialdata as spd
import sopa

from argparse import ArgumentParser as AP
from spatialhub.scripts.utils import die, print_arguments, RESET, GREEN

sopa.settings.auto_save_on_disk = False
sopa.settings.parallelization_backend = "dask"
sopa.settings.dask_client_kwargs["n_workers"] = 2

PIXEL_SIZE = 0.120280945  # µm (CosMx default)
DP = 15 / PIXEL_SIZE  # Diameter in pixels, assuming 15 µm cell diameter


def main():
    p = AP(description="Add CellPose segmentation shapes to the Zarr object.")
    p.add_argument("--zarr-in", required=True, help="Input Zarr directory.")
    p.add_argument("--ch-nuc", required=True, help="Channel name for nucleus.")
    p.add_argument("--ch-cyto", default=None, help="Channel name for cytoplasm.")
    p.add_argument("--mdl-path", required=True, help="Path to CellPose model.")
    p.add_argument("--dp", default=DP, type=float, help="Prior cell diameter (pixels).")
    p.add_argument("--gpu", action="store_true", help="Use GPU for CellPose.")
    p.add_argument("--img-key", default="image", help="Zarr key for the image data.")
    p.add_argument("--flow-thr", default=0.4, type=float, help="Flow threshold.")
    p.add_argument("--pb-thr", default=0, type=float, help="Cellprob threshold.")
    p.add_argument("--clip-limit", default=0.2, type=float, help="CLAHE clip limit.")
    p.add_argument("--gaussian-sigma", default=1, type=float, help="Gaussian sigma.")
    args = p.parse_args()
    print_arguments(args)

    sdata = spd.read_zarr(args.zarr_in)
    if args.img_key not in sdata:
        raise KeyError(f"zarr key '{args.img_key}' not found in dataset")
    if args.ch_nuc not in sopa.utils.get_channel_names(sdata):
        raise KeyError(f"channel '{args.ch_nuc}' not found in Zarr dataset")
    if args.ch_cyto and args.ch_cyto not in sopa.utils.get_channel_names(sdata):
        raise KeyError(f"channel '{args.ch_cyto}' not found in Zarr dataset")
    if "cellpose_segments" in sdata:
        raise KeyError(f"zarr key 'cellpose_segments' already exists in dataset")

    if args.ch_cyto:
        channels = [args.ch_cyto, args.ch_nuc]
    else:
        channels = args.ch_nuc
    print(f"Using channel {channels} for CellPose segmentation.")

    dim = sdata[args.img_key].shape[1:3]
    sopa.make_image_patches(sdata, patch_width=max(dim), patch_overlap=0)
    sopa.segmentation.cellpose(
        sdata,
        channels=channels,
        diameter=args.dp,
        pretrained_model=args.mdl_path,
        gpu=args.gpu,
        image_key=args.img_key,
        flow_threshold=args.flow_thr,
        cellprob_threshold=args.pb_thr,
        clip_limit=args.clip_limit,
        gaussian_sigma=args.gaussian_sigma,
        key_added="cellpose_segments",
    )
    sdata.write_element("cellpose_segments")

    print(f"{GREEN}Successfully added CellPose segmentation shapes.{RESET}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        die("Failed to add segmentation shapes.")
