from pathlib import Path

import anndata as ad


def load_matching_h5ads(
    h5ad_dir: Path,
    sample_id: str,
    points_from: str,
    shapes: str,
    coords: str,
    agg_func: str,
) -> ad.AnnData:
    """Load all H5ADs matching a set of required elements for a given sample."""

    if not h5ad_dir.exists():
        raise FileNotFoundError(f"H5AD directory not found: {h5ad_dir}")
    if not h5ad_dir.is_dir():
        raise NotADirectoryError(f"Expected an H5AD directory: {h5ad_dir}")

    # It only makes sense to concatenate H5ADs derived from the same segmentation mask (shapes),
    # the same coords system and using the same aggregation function
    # => let's pre-filter the universer of H5ADs to only those that match these criteria
    matching_h5ads = sorted(
        path 
        for path in h5ad_dir.glob("*.h5ad")
        if shapes in path.stem and coords in path.stem and agg_func in path.stem
    )

    if len(matching_h5ads) == 0:
        raise FileNotFoundError(
            f"No H5AD files found for the given segmentation mask (shapes), coordinates system and aggregation function."
        )

    # Parse comma-separated list from YAML
    element_list = [e.strip() for e in points_from.split(",") if e.strip()]

    h5ad_files = []

    for element in element_list:

        matches = sorted(
            path
            for path in matching_h5ads
            if sample_id in path.stem and element in path.stem
        )

        if len(matches) == 0:
            raise FileNotFoundError(
                f"No H5AD found for sample_id='{sample_id}', "
                f"element='{element}'."
            )

        if len(matches) > 1:
            raise ValueError(
                f"Multiple H5ADs found for sample_id='{sample_id}', "
                f"element='{element}': {[m.name for m in matches]}"
            )

        h5ad_files.append(matches[0])

    print(f"Found {len(h5ad_files)} H5AD file(s):")
    for path in h5ad_files:
        print(f"  {path.name}")

    # Loading adatas as a list and validating that they all have the same sample_id
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

        adatas.append(current)

    return adatas


def concatenate_adatas(adatas: list[ad.AnnData]) -> ad.AnnData:
    """Concatenate a list of AnnData objects along the var (probe) axis."""

    adata = ad.concat(
        adatas,
        axis="var",
        join="outer",
        merge="unique",
        uns_merge="unique",
    )

    print(
        f"Concatenated AnnData shape: "
        f"{adata.n_obs} cells x {adata.n_vars} probes"
    )

    return adata


def main():
    p = AP(description="Concatenate AnnData objects across the *probe* axis for one sample,"\n,
                        "segmentation mask (shapes), coordinates system and aggregation function.")
    p.add_argument("h5adout", help="Path to output h5ad file.")
    p.add_argument("--h5ad_dir", help="Path to input h5ad files.")
    p.add_argument("--sample_id", help="Sample ID for which to concatenate H5AD files.")
    p.add_argument("--points-from", required=True, help="List of Points to be aggregated.")
    p.add_argument("--shapes", required=True, help="Shapes from which the original H5AD was derived.")
    p.add_argument(
        "--agg-func",
        default="sum",
        choices=["mean", "sum", "median"],
        help="Aggregation function used for original aggregation.",
    )
    p.add_argument(
        "--coords",
        default="global",
        help="Target coordinate system for the original aggregated AnnData.",
    )
    args = p.parse_args()
    print_arguments(args)

    adatas = load_matching_h5ads(
        h5ad_dir=args.h5ad_dir,
        sample_id=args.sample_id,
        points_from=args.points_from,
        shapes=args.shapes,
        coords=args.coords,
        agg_func=args.agg_func,
    )

    adata = concatenate_adatas(adatas)
    adata.write(args.h5adout, compression="gzip")

    print(f"{GREEN}Successfully concatenated AnnDatas.{RESET}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        die("Failed to concatenate AnnData.")
