# Util to initialise project-level observation metadata
# from SpatialData Shapes (observations, whether segments or bins)

from argparse import ArgumentParser
from pathlib import Path

import geopandas as gpd
import pandas as pd
import spatialdata as sd
import sopa


def transform_shapes(sdata, shape_name, coords):
    """Return a Shapes GeoDataFrame expressed in the requested coordinate system."""
    transformed = sdata.transform_element_to_coordinate_system(
        element_name=shape_name, target_coordinate_system=coords
    )
    # SpatialData releases may return the transformed element or a SpatialData object.
    shapes = transformed[shape_name] if isinstance(transformed, sd.SpatialData) else transformed
    if not isinstance(shapes, gpd.GeoDataFrame):
        raise TypeError(f"Expected transformed Shapes GeoDataFrame for {shape_name!r}")
    return shapes


def init_sample_metadata(path, shapes_name, fov_name, coords):
    sdata = sd.read_zarr(path)
    shapes = transform_shapes(sdata, shapes_name, coords)
    if not shapes.index.is_unique:
        raise ValueError(f"Non-unique shape identifiers in {path}")

    # Extract minimal centroids coordinates and area for each shape.
    centroids = shapes.geometry.centroid
    meta = shapes.drop(columns="geometry").copy()
    meta.index = meta.index.map(str)
    meta.index.name = "instance_id"
    meta = meta.reset_index()
    
    if "sample_id" in meta.columns:
        meta = meta.drop(columns="sample_id")
    meta.insert(0, "sample_id", path.stem)
    
    if any(c in meta.columns for c in ("array_col", "array_row", "FOV", "area")):
        raise ValueError(f"Reserved metadata columns already exist in {shapes_name!r} for {path}")
    
    meta["array_col"] = centroids.x.to_numpy()
    meta["array_row"] = centroids.y.to_numpy()
    meta["area"] = shapes.geometry.area.to_numpy()

    # Assign each segment to its original FOV, using pseudoSNR to resolve overlaps
    # WARNING: will need adaptation for Xenium
    if fov_name is None:
        return meta

    fovs = transform_shapes(sdata, fov_name, coords)
    if not fovs.index.is_unique:
        raise ValueError(f"Non-unique FOV identifiers in {path}")
    if "pseudoSNR" not in fovs.columns or fovs["pseudoSNR"].isna().any():
        raise ValueError(f"Missing FOV pseudoSNR values in {path}")

    instances = gpd.GeoDataFrame(
        {"_row": range(len(shapes))}, geometry=centroids.to_numpy(), crs=shapes.crs
    )
    fov_bounds = fovs[["geometry", "pseudoSNR"]].copy()
    fov_bounds.index.name = None
    
    hits = gpd.sjoin(instances, fov_bounds, how="left", predicate="intersects")
    hits = hits.dropna(subset=["index_right"])
    
    # Highest pseudoSNR wins; stable string tie-break for equal scores.
    hits["_fov"] = hits["index_right"].astype(str)
    hits = hits.sort_values(
        ["_row", "pseudoSNR", "_fov"], ascending=[True, False, True], kind="stable"
    ).drop_duplicates("_row")
    
    meta["FOV"] = pd.Series(hits["index_right"].to_numpy(), index=hits["_row"].to_numpy()).reindex(range(len(meta))).to_numpy()
    if meta["FOV"].isna().any():
        raise ValueError(f"{meta['FOV'].isna().sum()} segments have no FOV match in {path}")
    
    return meta


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--zarr-dir", required=True)
    parser.add_argument("--shapes-by", required=True)
    parser.add_argument("--fov-shapes", default=None,
                        help="Optional FOV Shapes element for pseudoSNR-based assignment")
    parser.add_argument("--coords", default="global")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    out = Path(args.out)
    if out.exists():
        raise FileExistsError(f"Metadata already exists: {out}")
    paths = sorted(p for p in Path(args.zarr_dir).glob("*.zarr") if p.is_dir())
    if not paths:
        raise ValueError(f"No sample Zarr directories found in {args.zarr_dir}")
    tables = []
    for path in paths:
        print(f"Processing {path.name}", flush=True)
        tables.append(init_sample_metadata(path, args.shapes_by, args.fov_shapes, args.coords))
    result = pd.concat(tables, ignore_index=True)
    if result.duplicated(["sample_id", "instance_id"]).any():
        raise ValueError("Duplicate (sample_id, instance_id) pairs")
    out.parent.mkdir(parents=True, exist_ok=True)
    # Exclusive creation prevents accidental overwriting of central annotations.
    with out.open("xb") as handle:
        result.to_parquet(handle, index=False)
    print(f"Wrote {len(result):,} instances from {len(paths)} samples to {out}")


if __name__ == "__main__":
    main()
