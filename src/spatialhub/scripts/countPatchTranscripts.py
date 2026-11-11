
#!/usr/bin/env python3

"""Count transcripts per SOPA patch without repeated spatial queries."""

import argparse
from pathlib import Path

import dask
import numpy as np
import pandas as pd
import spatialdata as sd


def count_partition(df, bounds):
    """Count transcripts in all patch bounding boxes for one partition."""

    x = df["x"].to_numpy()
    y = df["y"].to_numpy()

    counts = np.zeros(len(bounds), dtype=np.int64)

    for i, (xmin, ymin, xmax, ymax) in enumerate(bounds):
        counts[i] = np.count_nonzero(
            (x >= xmin) & (x <= xmax) &
            (y >= ymin) & (y <= ymax)
        )

    return counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-zarr", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    sdata = sd.read_zarr(args.input_zarr)

    patches = sdata["transcripts_patches"]
    transcripts = sdata.points["transcripts"]

    if patches.empty:
        raise ValueError("No transcript patches found.")

    patch_indices = patches.index.tolist()

    # SOPA patch bounding boxes
    bounds = np.asarray(
        patches["bboxes"].tolist(),
        dtype=float,
    )

    # Process each Dask partition once.
    tasks = [
        dask.delayed(count_partition)(partition, bounds)
        for partition in transcripts.to_delayed()
    ]

    partition_counts = dask.compute(
        *tasks,
        scheduler="threads",
        num_workers=args.workers,
    )

    counts = np.sum(
        np.stack(partition_counts),
        axis=0,
    )

    result = pd.DataFrame({
        "patch_index": patch_indices,
        "n_transcripts": counts,
    })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False)

    print(result.to_string(index=False))
    print(f"\nTotal patches: {len(result)}")
    print(f"Minimum transcripts: {result.n_transcripts.min()}")
    print(f"Maximum transcripts: {result.n_transcripts.max()}")


if __name__ == "__main__":
    main()
