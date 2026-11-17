import os
import numpy as np
import spatialdata as spd
from spatialdata.transformations import (
    Affine,
    get_transformation,
    set_transformation,
)
from scipy import sparse
from spatialdata.models import PointsModel

from argparse import ArgumentParser as AP
from pathlib import Path
from spatialhub.scripts.utils import die, print_arguments, RESET, GREEN


def affine_matrix(transform, axes):
    """Return transform as an affine matrix on the requested axes."""
    return transform.to_affine_matrix(
        input_axes=axes,
        output_axes=axes,
    )


def validate_transform(sdata, points_from, shapes_by, coords):
    """
    Prepare points for aggregation against shapes.
    Preserve valid transformations, reducing the point
    transformation from 3D to 2D only when necessary.
    """
    
    shapes = sdata.shapes[shapes_by]
    points = sdata.points[points_from]

    # Check the gene/probe/feature name is correctly registered
    point_attrs = points.attrs.get(PointsModel.ATTRS_KEY, {})
    feature_key = point_attrs.get(PointsModel.FEATURE_KEY)

    if feature_key is None:
        raise ValueError(
            "No feature key registered for points element "
            f"'{args.points_from}'. Cannot generate a gene-by-cell matrix."
        )

    if feature_key not in points.columns:
        raise ValueError(
            f"Registered feature column '{feature_key}' "
            "is missing from the points table."
        )

    print(f"Using transcript feature column: {feature_key}")

    # Determine dimensionality from the columns in the Points (transcripts) table.
    point_axes = tuple(
        axis for axis in ("x", "y", "z")
        if axis in points.columns
    )
    if point_axes not in (("x", "y"), ("x", "y", "z")):
        raise ValueError(f"Unsupported point axes: {point_axes}")

    # Determine dimensionality from the geometries in the Shapes (cells) GeoPandas DataFrame.
    geometry = shapes.geometry

    if geometry.isna().any():
        raise ValueError("Shapes contain missing geometries.")

    if geometry.is_empty.any():
        raise ValueError("Shapes contain empty geometries.")

    #if not geometry.is_valid.all():
    #    raise ValueError("Shapes contain invalid geometries.")

    has_z = shapes.geometry.has_z
    if has_z.any() and not has_z.all():
        raise ValueError(
            "Shapes contain mixed 2D and 3D geometries."
        )
    shape_axes = ("x", "y", "z") if has_z.all() else ("x", "y")
    
    point_transform = get_transformation(points, to_coordinate_system=coords)
    shape_transform = get_transformation(shapes, to_coordinate_system=coords)

    # ------------------------------------------------------------------
    # Same geometric dimensionality - fine, moving on!
    # might need some refinement later (edge cases with axes incompatibilities)
    # ------------------------------------------------------------------
    if point_axes == shape_axes:
        return points
    
    # ------------------------------------------------------------------
    # 2D points + 3D shapes:
    # => don't silently discard the shape z dimension.
    # ------------------------------------------------------------------
    elif point_axes == ("x", "y") and shape_axes == ("x", "y", "z"):

        raise ValueError(
            "Cannot aggregate 2D points by 3D shapes "
            "without an explicit 3D->2D projection."
        )
    
    # ------------------------------------------------------------------
    # 3D points + 2D shapes:
    # compare the x/y component of the transformations.
    # ------------------------------------------------------------------
    elif point_axes == ("x", "y", "z") and shape_axes == ("x", "y"):

        matrix_3d = affine_matrix(
            point_transform,
            ("x", "y", "z"),
        )

        # The x/y output must not depend on z.
        #
        # Matrix layout:
        #
        #   x'   [ ... ... M[x,z] tx ] [x]
        #   y' = [ ... ... M[y,z] ty ] [y]
        #   z'   [ ... ... ...     tz ] [z]
        #
        # If M[x,z] or M[y,z] are non-zero, z cannot safely be discarded.
        xy_from_z = matrix_3d[:2, 2]

        if not np.allclose(xy_from_z, 0):
            raise ValueError(
                "Cannot reduce 3D point transformation to 2D: "
                "x/y coordinates depend on z. "
                f"x/y <- z terms: {xy_from_z}"
            )

        # Extract the 2D affine transformation.
        matrix_2d = matrix_3d[np.ix_(
            [0, 1, 3],
            [0, 1, 3],
        )]

        # Work on a copy; Keep the original SpatialData element unchanged.
        points = points.copy()
        points.attrs["transform"] = dict(points.attrs["transform"])  # 'detaches' the copy from its source

        set_transformation(
            points,
            Affine(
                matrix_2d,
                input_axes=("x", "y"),
                output_axes=("x", "y"),
            ),
            to_coordinate_system=coords,
        )

        print("Reduced 3D point transformation to 2D.")

        return points


def filter_points(points):
    # This is for Xenium data
    if "qv" in points.columns:
        points = points[points['qv'] > 20]  # add option to modify this filtering threshold
    return points


def main():
    p = AP(description="Aggregate to AnnData object from SpatialData Zarr file.")
    p.add_argument("h5adout", help="Output h5ad file path.")
    p.add_argument("--in-zarr", required=True, help="Input Zarr file path.")
    p.add_argument("--points-from", required=True, help="Points to be aggregated.")
    p.add_argument("--shapes-by", required=True, help="Shapes to aggregate by.")
    p.add_argument(
        "--agg-func",
        default="sum",
        choices=["mean", "sum", "median"],
        help="Aggregation function.",
    )
    p.add_argument(
        "--coords",
        default="global",
        help="Target coordinate system for the aggregated AnnData.",
    )
    args = p.parse_args()
    print_arguments(args)

    if os.path.exists(args.h5adout):
        raise FileExistsError(f"Output file {args.h5adout} already exists.")

    sdata = spd.read_zarr(args.in_zarr)
    points = validate_transform(
        sdata,
        points_from=args.points_from,
        shapes_by=args.shapes_by,
        coords=args.coords,
    )
    points = filter_points(points)

    aggregated = sdata.aggregate(
        values=points,
        by=args.shapes_by,
        agg_func=args.agg_func,
        target_coordinate_system=args.coords,
        deepcopy=False,
    )

    adata = aggregated.tables["table"]
    
    # Add spatial metadata to adata.obs
    shapes = aggregated.shapes[args.shapes_by]
    # Convert shape identifiers to strings to match AnnData obs_names.
    centroids = shapes.geometry.centroid

    if not centroids.index.is_unique:
        raise ValueError("Shape identifiers are not unique.")

    centroids = centroids.reindex(adata.obs_names)

    if centroids.isna().any():
        raise ValueError(
            "Some AnnData cells have no matching shape centroid."
        )
    
    adata.obs["array_col"] = centroids.x.to_numpy()
    adata.obs["array_row"] = centroids.y.to_numpy()

    adata.obsm["spatial"] = adata.obs[
        ["array_col", "array_row"]
    ].to_numpy()
    
    # Add other essential metadata to adata.obs
    sample_id = Path(args.in_zarr).stem
    adata.obs["sample_id"] = sample_id
    
    adata.layers["counts"] = adata.X.copy()  # Ensure counts layer is present
    
    # Store provenance of the aggrgated counts matrix in adata.uns
    # to avoid loss of information in case of file name changes.
    adata.uns["counts_mtx_source"] = {
        "input_zarr": str(args.in_zarr),
        "sample_id": sample_id,
        "points_from": [
            point.strip()
            for point in args.points_from.split(",")
            if point.strip()
        ],
        "shapes_by": args.shapes_by,
        "agg_func": args.agg_func,
        "coords_system": args.coords,
    }

    # Perform some sanity checks!
    X = adata.X
    values = X.data if sparse.issparse(X) else np.asarray(X).ravel()

    if not np.all(np.isfinite(values)):
        raise ValueError("Aggregated matrix contains non-finite values.")

    if np.any(values < 0):
        raise ValueError("Aggregated matrix contains negative counts.")

    if not np.all(values == np.floor(values)):
        raise ValueError("Aggregated matrix contains non-integer counts.")

    print(f"Aggregated matrix: {adata.n_obs} cells × {adata.n_vars} genes")

    # Write final AnnData H5AD file
    adata.write(args.h5adout, compression="gzip")

    print(f"{GREEN}Successfully aggregated AnnData.{RESET}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        import traceback
        traceback.print_exc()
        die("Failed to aggregate AnnData.")
