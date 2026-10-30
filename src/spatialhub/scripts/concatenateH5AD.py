import anndata as ad

from argparse import ArgumentParser as AP
from pathlib import Path
from spatialhub.scripts.utils import die, print_arguments, RESET, GREEN


def concatenate_adatas(adatas: list[ad.AnnData]) -> ad.AnnData:
    """
    Concatenate AnnData objects along the var (probe) axis.

    All counts_mtx_source metadata must be identical between inputs
    except for points_from, which is combined across inputs.
    """

    if len(adatas) < 2:
        raise ValueError("At least two input matrices are required for concatenation.")

    # Check whether AnnDatas to concatenate wwere generated using spatialhub or not
    # and save provenance accordingly if applicable
    sources = [a.uns.get("counts_mtx_source") for a in adatas]
    if any(s is not None for s in sources) and not all(s is not None for s in sources):
        raise ValueError("Inconsistent availability of 'counts_mtx_source' uns metadata")
    spatialhub = all(s is not None for s in sources)

    if spatialhub:
        provenance = {k: v for k, v in sources[0].items() if k != "points_from"}
        for source in sources[1:]:
            if {k: v for k, v in source.items() if k != "points_from"} != provenance:
                raise ValueError("Incompatible 'counts_mtx_source' metadata")

    # Initiate list to build on and object to check against
    first = adatas[0]
    observed_features = set()
    aligned_adatas = []
    
    for i, current in enumerate(adatas, start=1):
        if not current.obs_names.is_unique or not current.var_names.is_unique:
            raise ValueError(
                f"AnnData {i} has duplicate observation or feature identifiers"
            )
        if not current.obs_names.equals(first.obs_names):
            if set(current.obs_names) != set(first.obs_names):
                raise ValueError(
                    f"AnnData {i} does not have the same observations as AnnData 1"
                )
            current = current[first.obs_names, :].copy()

        duplicates = observed_features.intersection(current.var_names)
        if duplicates:
            raise ValueError(
                f"Duplicate feature identifiers across matrices: "
                f"{sorted(duplicates)[:5]}"
            )
        observed_features.update(current.var_names)
        
        aligned_adatas.append(current)

    # Concatenate aligned AnnDatas
    result = ad.concat(
        aligned_adatas, axis=1, join="inner", merge="first", uns_merge="first"
    )

    # Retrieve minimal metadata
    if spatialhub:
        required_obs = ["sample_id", "array_col", "array_row"]
        missing = set(required_obs) - set(first.obs.columns)
        if missing:
            raise ValueError(
                f"First AnnData is missing required spatialhub metadata: {sorted(missing)}"
            )

        result.obs = first.obs[required_obs].copy()
        result.uns["counts_mtx_source"] = {
            **provenance,
            "points_from": sum((s["points_from"] for s in sources), []),
        }
    else:
        result.obs = first.obs.copy()
        print(
            "No spatialhub provenance detected. "
            "Retaining only observation metadata from the first AnnData.\n"
            f"Available columns: {list(first.obs.columns)}"
        )

    return result


def main():
    p = AP(
        description=(
            "Concatenate AnnData objects across the probe axis for one sample, "
            "segmentation mask, coordinate system and aggregation function."
        )
    )
    p.add_argument("h5adout", help="Path to output H5AD file.")
    p.add_argument("inputs", nargs="+")
    args = p.parse_args()
    print_arguments(args)

    output = Path(args.h5adout)
    if output.exists():
        raise FileExistsError(output)

    # Since the first AnnData provides the authoritative metadata, the order matters.
    inputs = sorted(args.inputs)
    print(f"\n{GREEN}AnnData input order (first supplies .obs metadata):{RESET}")
    for i, path in enumerate(inputs, start=1):
        print(f"  {i}. {path}")
    
    result = concatenate_adatas([ad.read_h5ad(path) for path in inputs])

    output.parent.mkdir(parents=True, exist_ok=True)
    result.write_h5ad(output, compression="gzip")

    print(f"{GREEN}Successfully concatenated AnnDatas.{RESET}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        die("Failed to concatenate AnnDatas.")
