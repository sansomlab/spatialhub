import argparse
import spatialdata as spd
import spatialdata_plot
import matplotlib.pyplot as plt

if __name__ == "__main__":
    print("Parsing arguments")
    parser = argparse.ArgumentParser(
        description="Plot images of a Zarr object to check."
    )
    parser.add_argument(
        "--zarr", required=True, type=str, help="Path to the input Zarr file."
    )
    parser.add_argument(
        "--output", required=True, type=str, help="Path to the output plot image."
    )
    args = parser.parse_args()

    sdata = spd.read_zarr(args.zarr)

    axes = plt.subplots(1, 4, figsize=(20, 5))[1]

    sdata.pl.render_images("image", scale="scale2").pl.show(ax=axes[0])
    axes[0].invert_yaxis()
    axes[0].set_title("image")

    sdata.pl.render_images("image", scale="scale2").pl.render_shapes(
        "shapes",
        outline_color="yellow",
        outline_width=1,
        outline_alpha=1,
        fill_alpha=0,
    ).pl.show(ax=axes[1])
    axes[1].invert_yaxis()
    axes[1].set_title("segmentations")

    sdata.pl.render_images("image", scale="scale2").pl.render_points(
        "tx_main", color="yellow", size=0.05, alpha=0.2
    ).pl.show(ax=axes[2])
    axes[2].invert_yaxis()
    axes[2].set_title("tx_main")

    sdata.pl.render_images("image", scale="scale2").pl.render_points(
        "tx_ctrl", color="yellow", size=0.05, alpha=0.6
    ).pl.show(ax=axes[3])
    axes[3].invert_yaxis()
    axes[3].set_title("tx_ctrl")

    plt.tight_layout()
    plt.savefig(args.output, bbox_inches="tight", dpi=300)
