import anndata as ad

from argparse import ArgumentParser as AP
from pathlib import Path
from spatialhub.scripts.utils import die, print_arguments, RESET, GREEN


def load_matching_h5ads(
    h5ad_dir: Path,
    sample_id: str,
    points_from: str,
    shapes: str,
    coords: str,
    agg_func: str,
) -> list[ad.AnnData]:
    """Load all H5ADs matching the specified aggregation configuration for a sample."""

    h5ad_dir = Path(h5ad_dir)

    if not h5ad_dir.exists():
        raise FileNotFoundError(f"H5AD directory not found: {h5ad_dir}")
    if not h5ad_dir.is_dir():
        raise NotADirectoryError(f"Expected an H5AD directory: {h5ad_dir}")

    elements = [e.strip() for e in points_from.split(",") if e.strip()]

    if not elements:
        raise ValueError("No Points elements specified in --points-from.")

    h5ad_files = [
        h5ad_dir
        / f"{element}.{shapes}.{agg_func}.{coords}"
        / f"{sample_id}.h5ad"
        for element in elements
    ]

    missing = [path for path in h5ad_files if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Required H5AD file(s) not found:\n"
            + "\n".join(str(path) for path in missing)
        )

    print(f"Found {len(h5ad_files)} H5AD file(s):")
    for path in h5ad_files:
        print(f"  {path}")

    adatas = []

    for path in h5ad_files:
        current = ad.read_h5ad(path)

        if "sample_id" not in current.obs.columns:
            raise ValueError(
                f"Required AnnData.obs column 'sample_id' not found in {path}"
            )

        sample_values = current.obs["sample_id"].dropna().astype(str).unique()

        if len(sample_values) != 1 or sample_values[0] != sample_id:
            raise ValueError(
                f"H5AD {path} does not contain exactly one sample_id "
                f"matching '{sample_id}'. Found: {sample_values.tolist()}"
            )

        if "counts_mtx_source" not in current.uns:
            raise ValueError(
                f"Required adata.uns['counts_mtx_source'] not found in {path}"
            )

        adatas.append(current)

    return adatas


def concatenate_adatas(adatas: list[ad.AnnData]) -> ad.AnnData:
    """
    Concatenate AnnData objects along the var (probe) axis.

    All counts_mtx_source metadata must be identical between inputs
    except for points_from, which is combined across inputs.
    """

    if not adatas:
        raise ValueError("No AnnData objects provided for concatenation.")

    sources = [
        dict(adata.uns["counts_mtx_source"])
        for adata in adatas
    ]

    # All provenance fields except points_from must be identical.
    common_source = {
        key: value
        for key, value in sources[0].items()
        if key != "points_from"
    }

    for i, source in enumerate(sources[1:], start=2):
        current = {
            key: value
            for key, value in source.items()
            if key != "points_from"
        }

        if current != common_source:
            raise ValueError(
                "Input H5ADs have inconsistent counts_mtx_source metadata "
                f"(input 1 vs input {i}):\n"
                f"Input 1: {common_source}\n"
                f"Input {i}: {current}"
            )

    # Combine points_from values, preserving input order.
    points_from = []

    for source in sources:
        if "points_from" not in source:
            raise ValueError(
                "Required 'points_from' field missing from "
                "adata.uns['counts_mtx_source']."
            )

        points = source["points_from"]

        if isinstance(points, str):
            points = [points]

        for point in points:
            point = str(point)
            if point not in points_from:
                points_from.append(point)

    adata = ad.concat(
        adatas,
        axis="var",
        join="outer",
        merge="unique",
        uns_merge="unique",
    )

    # Replace the differing input provenance with the validated,
    # combined provenance for the new count matrix.
    adata.uns["counts_mtx_source"] = {
        **common_source,
        "points_from": points_from,
    }

    print(
        f"Concatenated AnnData shape: "
        f"{adata.n_obs} cells x {adata.n_vars} probes"
    )
    print("Count matrix source:")
    print(adata.uns["counts_mtx_source"])

    return adata


def main():
    p = AP(
        description=(
            "Concatenate AnnData objects across the probe axis for one sample, "
            "segmentation mask, coordinate system and aggregation function."
        )
    )
    p.add_argument("h5adout", help="Path to output H5AD file.")
    p.add_argument(
        "--h5ad_dir",
        required=True,
        help="Path to input H5AD files.",
    )
    p.add_argument(
        "--sample_id",
        required=True,
        help="Sample ID for which to concatenate H5AD files.",
    )
    p.add_argument(
        "--points-from",
        required=True,
        help="Comma-separated list of Points elements to concatenate.",
    )
    p.add_argument(
        "--shapes",
        required=True,
        help="Shapes from which the original H5ADs were derived.",
    )
    p.add_argument(
        "--agg-func",
        default="sum",
        choices=["mean", "sum", "median"],
        help="Aggregation function used for original aggregation.",
    )
    p.add_argument(
        "--coords",
        default="global",
        help="Coordinate system used for original aggregation.",
    )

    args = p.parse_args()
    print_arguments(args)

    h5adout = Path(args.h5adout)

    if h5adout.exists():
        raise FileExistsError(f"Output file {h5adout} already exists.")

    adatas = load_matching_h5ads(
        h5ad_dir=args.h5ad_dir,
        sample_id=args.sample_id,
        points_from=args.points_from,
        shapes=args.shapes,
        coords=args.coords,
        agg_func=args.agg_func,
    )

    adata = concatenate_adatas(adatas)
    adata.write(h5adout, compression="gzip")

    print(f"{GREEN}Successfully concatenated AnnDatas.{RESET}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        die("Failed to concatenate AnnData.")
