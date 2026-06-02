import os
import argparse
import pandas as pd


def parse_grid_smpgrp(fovpos_smp, width, height):
    xmin_px = fovpos_smp["x_global_px"].min()
    ymax_px = fovpos_smp["y_global_px"].max()

    out = pd.DataFrame(index=fovpos_smp.index)
    out["FOV"] = fovpos_smp.index
    out["x_igrid"] = ((fovpos_smp["x_global_px"] - xmin_px) / width + 0.5).astype(int)
    out["y_igrid"] = ((ymax_px - fovpos_smp["y_global_px"]) / height + 0.5).astype(int)

    return out.reset_index(drop=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--fovpos-csv", required=True, help="CSV for FOV positions.")
    p.add_argument("--raw-meta", required=True, help="TSV for raw metadata.")
    p.add_argument("--slide", required=True, help="Only process this slide ID.")
    p.add_argument("--width", type=int, default=4256, help="FOV width.")
    p.add_argument("--height", type=int, default=4256, help="FOV height.")
    p.add_argument("--tsv-out", default=None, help="Write to TSV; defaults to stdout.")
    args = p.parse_args()

    meta = pd.read_csv(args.raw_meta, sep="\t")
    assert "sample_id" in meta.columns, "metadata needs a 'sample_id' column."
    assert "slide_id" in meta.columns, "metadata needs a 'slide_id' column."
    assert "fov_sequence" in meta.columns, "metadata needs a 'fov_sequence' column."
    meta = meta.loc[meta["slide_id"] == args.slide].copy()
    meta.set_index("sample_id", inplace=True)

    fovpos = pd.read_csv(args.fovpos_csv)
    if "FOV" in fovpos.columns:
        fovpos.set_index("FOV", inplace=True)
    elif "fov" in fovpos.columns:
        fovpos.set_index("fov", inplace=True)
    else:
        raise ValueError("fovpos CSV needs a column named either 'FOV' or 'fov'.")
    assert "x_global_px" in fovpos.columns, "fovpos CSV needs 'x_global_px' column."
    assert "y_global_px" in fovpos.columns, "fovpos CSV needs 'y_global_px' column."
    fovpos = fovpos[["x_global_px", "y_global_px"]]

    # make dict for mapping FOV to sample ID, for assigning sample ID to fovpos
    fov2smp = (
        meta["fov_sequence"]
        .str.split(",")
        .explode()
        .reset_index()
        .astype({"fov_sequence": int})
        .set_index("fov_sequence")["sample_id"]
        .to_dict()
    )
    fovpos["sample_id"] = fovpos.index.map(fov2smp)

    # make FOV grid coordinates for each sample, and assign to fovpos
    fovgrids = (
        fovpos.groupby("sample_id")[["x_global_px", "y_global_px"]]
        .apply(lambda df: parse_grid_smpgrp(df, args.width, args.height))
        .set_index("FOV")
    )
    fovpos = (
        pd.concat([fovpos, fovgrids], axis=1, join="inner", verify_integrity=True)
        .reset_index()
        .rename(columns={"index": "FOV"})
    )

    # make dicts for mapping sample ID to FOV sequence, width, and height, for assigning to metadata
    smp2width = (fovpos.groupby("sample_id")["x_igrid"].max() + 1).to_dict()
    smp2height = (fovpos.groupby("sample_id")["y_igrid"].max() + 1).to_dict()
    smp2fovseq = {}
    for smp in fovpos["sample_id"].unique():
        fovlst = []
        smpgrid2fov = (
            fovpos.loc[fovpos["sample_id"] == smp]
            .set_index(["x_igrid", "y_igrid"])["FOV"]
            .to_dict()
        )
        for y_igrid in range(smp2height[smp]):
            for x_igrid in range(smp2width[smp]):
                fovlst += [str(smpgrid2fov.get((x_igrid, y_igrid), "blank"))]
        assert smp not in smp2fovseq
        smp2fovseq[smp] = ",".join(fovlst)

    # assign FOV sequence, width, and height to metadata
    meta["fov_sequence"] = meta.index.map(smp2fovseq)
    meta["fov_width"] = meta.index.map(smp2width)
    meta["fov_height"] = meta.index.map(smp2height)

    if args.tsv_out is not None:
        meta.reset_index().to_csv(args.tsv_out, sep="\t", index=False)
    else:
        print(meta.reset_index().to_csv(sep="\t", index=False))
