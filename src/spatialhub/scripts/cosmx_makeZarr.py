import os
import numpy as np
import pandas as pd
import geopandas as gpd
import tifffile
import spatialdata as spd
from scipy.spatial import cKDTree

from argparse import ArgumentParser as AP
from shapely import Polygon
from spatialhub.scripts.utils import die, print_arguments, RESET, YELLOW, GREEN


def transform_coordinates(points_df, fov2pos):
    points_df["x"] = points_df.apply(
        lambda r: fov2pos[r["FOV"]]["Position_X"]
        + r["x_global_px"]
        - fov2pos[r["FOV"]]["x_global_px"],
        axis=1,
    )
    points_df["y"] = points_df.apply(
        lambda r: fov2pos[r["FOV"]]["Position_Y"]
        + fov2pos[r["FOV"]]["y_global_px"]
        - r["y_global_px"],
        axis=1,
    )
    return points_df


# Quality metric to rank FOVs by pseudo-SNR (signal-to-noise ratio)
# used to filter transcripts in overlapping FOV regions
def calculate_pseudoSNR(points_df, rna_probes, neg_probes, alpha = 3):
    rna_probes = set(rna_probes)
    neg_probes = set(neg_probes)

    if not rna_probes:
        raise ValueError("RNA probe set must be provided")
    if rna_probes & neg_probes:
        raise ValueError("RNA and negative probe sets overlap")

    # Count RNA and negative-control transcripts per FOV
    counts = (
        points_df.assign(
            is_rna=points_df["target"].isin(rna_probes),
            is_neg=points_df["target"].isin(neg_probes),
        )
        .groupby("FOV", sort=False)[["is_rna", "is_neg"]]
        .sum()
        .rename(columns={
            "is_rna": "nCount_RNA_FOVsum",
            "is_neg": "nCount_neg_FOVsum",
        })
    )

    # Normalise by the number of distinct probes in the panel
    # + pseudocount alpha, to avoid infinite values on FOVs with no negative probe detection
    counts["FOV_SNR"] = counts["nCount_RNA_FOVsum"] / (counts["nCount_neg_FOVsum"] + alpha)

    return counts.reset_index()


def filter_overlap_tx(
    transcripts,
    fovpos,
    tile_width,
    tile_height,
):
    """
    Remove transcripts in overlapping FOV regions when another
    overlapping FOV has a higher SNR.

    Parameters
    ----------
    transcripts : pd.DataFrame
        Must contain FOV, x, y
    fovpos : pd.DataFrame
        Must contain FOV, Position_X, Position_Y, FOV_SNR.

    tile_width, tile_height : int
        Dimensions of each FOV image in pixels.

    Returns
    -------
    pd.DataFrame
        Filtered transcripts, preserving original columns and index.

    Notes
    -----
    - Uses ASHLAR's final stitched FOV positions.
    - Only examines transcripts in FOVs with overlapping neighbours.
    - Highest-SNR FOV owns each overlapping region.
    - SNR ties are resolved by retaining the lowest FOV identifier.
    - A transcript is removed if any higher-priority FOV covers
      its position, even if that FOV did not detect the molecule.
    """

    # Validate inputs
    fovpos = fovpos.reset_index(drop=False)  # Ensure FOV is a column, not index
    required_t = {"FOV", "x", "y"}
    required_f = {"FOV", "Position_X", "Position_Y", "FOV_SNR"}

    if missing := required_t - set(transcripts.columns):
        raise ValueError(f"Missing transcript columns: {missing}")
    if missing := required_f - set(fovpos.columns):
        raise ValueError(f"Missing FOV columns: {missing}")

    if not fovpos["FOV"].is_unique:
        raise ValueError("Duplicate FOV identifiers")

    if not (tile_width > 0 and tile_height > 0):
        raise ValueError("Tile dimensions must be positive")

    if not np.isfinite(
        fovpos[["Position_X", "Position_Y", "FOV_SNR"]]
        .to_numpy(dtype=float)
    ).all():
        raise ValueError("Invalid FOV positions or SNR values")

    if not np.isfinite(
        transcripts[["x", "y"]]
        .to_numpy(dtype=float)
    ).all():
        raise ValueError("Invalid transcript coordinates")

    if not transcripts["FOV"].isin(fovpos["FOV"]).all():
        raise ValueError("Transcripts contain unknown FOV identifiers")

    # Sort FOVs by priority: highest SNR first, then lowest FOV ID
    fovpos = fovpos.sort_values(
        ["FOV_SNR", "FOV"],
        ascending=[False, True],
    ).reset_index(drop=True)

    fovs = fovpos["FOV"].to_numpy()
    origins = fovpos[
        ["Position_X", "Position_Y"]
    ].to_numpy(dtype=float)

    # Identify potentially overlapping FOVs using a KD-tree
    tree = cKDTree(origins)

    # Maximum possible distance between origins of overlapping tiles
    radius = np.hypot(tile_width, tile_height)
    pairs = tree.query_pairs(radius, output_type="ndarray")

    # Group transcript row positions by FOV once
    groups = transcripts.groupby("FOV", sort=False).indices

    coords = transcripts[
        ["x", "y"]
    ].to_numpy(dtype=float)

    keep = np.ones(len(transcripts), dtype=bool)
    n_overlaps = 0

    for i, j in pairs:

        xi, yi = origins[i]
        xj, yj = origins[j]

        # Calculate the actual overlap rectangle
        xmin = max(xi, xj)
        xmax = min(xi + tile_width, xj + tile_width)
        ymin = max(yi, yj)
        ymax = min(yi + tile_height, yj + tile_height)

        if xmin >= xmax or ymin >= ymax:
            continue

        n_overlaps += 1

        # FOVs are sorted by descending SNR,
        # so the higher index has lower priority.
        loser = fovs[max(i, j)]

        rows = groups.get(loser)
        if rows is None or len(rows) == 0:
            continue

        # Only check transcripts not already removed
        rows = rows[keep[rows]]
        xy = coords[rows]

        inside = (
            (xy[:, 0] >= xmin)
            & (xy[:, 0] < xmax)
            & (xy[:, 1] >= ymin)
            & (xy[:, 1] < ymax)
        )

        keep[rows[inside]] = False

    print(f"Overlapping FOV pairs: {n_overlaps}")
    print(f"Input transcripts: {len(transcripts):,}")
    print(f"Removed transcripts: {(~keep).sum():,}")
    print(f"Retained transcripts: {keep.sum():,}")

    return transcripts.iloc[np.flatnonzero(keep)].copy()


def main():
    p = AP(description="Prepare a Zarr object from CosMx data.")
    # Required arguments
    p.add_argument("output", help="Output Zarr directory.")
    p.add_argument("--fov-csv", required=True, help="CosMx FOV position CSV file.")
    p.add_argument("--tile-size", required=True, help="Size of FOV tile.")
    p.add_argument("--tx-csv", required=True, help="CosMx transcript CSV file.")
    p.add_argument("--morph-tiff", required=True, help="Assembled/stitched TIFF image.")
    # Optional arguments
    p.add_argument("--pg-csv", default=None, help="CosMx polygon CSV file.")
    p.add_argument("--channels", default=None, help="Comma-separated channel names.")
    p.add_argument("--ctrl-regex", default='Negative|SystemControl', help="Regex for control probes, including positive and negative ones.")
    p.add_argument("--other-regex", default=None, help="Regex for other probes.")
    p.add_argument("--scales", default=None, help="Comma-separated image scales.")
    p.add_argument("--alpha", type=float, default=3.0, help="Pseudocount for SNR calculation (default: 3).")
    args = p.parse_args()
    print_arguments(args)

    if os.path.exists(args.output):
        raise FileExistsError(f"output file {args.output} already exists")
    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    print(f"Output Zarr directory will be created at '{args.output}'.")

    if not os.path.exists(args.fov_csv):
        raise FileNotFoundError(f"FOV position file '{args.fov_csv}' not found")
    if not os.path.exists(args.tx_csv):
        raise FileNotFoundError(f"Transcript file '{args.tx_csv}' not found")
    if not os.path.exists(args.morph_tiff):
        raise FileNotFoundError(f"TIFF file '{args.morph_tiff}' not found")
    
    if not np.isfinite(args.alpha) or args.alpha <= 0:
        raise ValueError("--alpha must be finite and strictly positive")

    
    ### Prepare Images ###

    img = tifffile.imread(args.morph_tiff)
    if args.channels:
        channels = [ch.strip() for ch in args.channels.split(",")]
    else:
        channels = list(range(img.shape[0]))

    sf = list(map(int, args.scales.split(","))) if args.scales else None
    print(f"Image shape: {img.shape}\nChannels: {channels}\nScale factors: {sf}")

    # Create a Image2DModel object for the TIFF image, with the specified channels and scale factors
    images = {
        "image": spd.models.Image2DModel.parse(img, c_coords=channels, scale_factors=sf)
    }


    ### Prepare Points ###

    # Read the FOV position CSV
    fovpos = pd.read_csv(args.fov_csv, index_col="FOV", header=0)
    missing = {"Position_X", "Position_Y"} - set(fovpos.columns)
    if missing:
        raise ValueError(f"missing {', '.join(sorted(missing))} in FOV table")
    fov2pos = fovpos[
        ["Position_X", "Position_Y", "x_global_px", "y_global_px"]
    ].to_dict(orient="index")

    # Read the transcript CSV and compute global coordinates for each transcript
    tx = pd.read_csv(args.tx_csv, header=0)
    if "fov" in tx.columns:
        tx.rename(columns={"fov": "FOV"}, inplace=True)
    missing = {"FOV", "x_global_px", "y_global_px", "target"} - set(tx.columns)
    if missing:
        raise ValueError(f"missing {', '.join(sorted(missing))} in tx table")
    tx = tx.loc[tx["FOV"].isin(fov2pos.keys())].copy()
    tx = transform_coordinates(tx, fov2pos)

    # Calculate SNR for each FOV based on RNA and negative control probes
    neg_probes = tx.loc[tx["target"].str.contains('Negative', na=False), "target"].unique().tolist()
    #print(f"Negative control probes: {set(neg_probes)}")
    
    patterns = [p for p in [args.ctrl_regex, args.other_regex] if p]
    exclude_regex = "|".join(f"(?:{p})" for p in patterns)

    rna_probes = tx.loc[
        (~tx["target"].str.contains(exclude_regex, na=False) if exclude_regex else pd.Series(True, index=tx.index)
        ) & ~tx["target"].isin(neg_probes),
        "target"
    ].unique().tolist()
    #print(f"RNA probes: {set(rna_probes)}")

    fov_snr = calculate_pseudoSNR(
        tx,
        rna_probes=rna_probes,
        neg_probes=neg_probes,
        alpha=args.alpha
    )
    # Append to FOV position file
    fovpos = fovpos.join(
         fov_snr.set_index("FOV")[["FOV_SNR"]].reindex(fovpos.index, fill_value=0.0),
         validate="one_to_one",
     )
    if fovpos["FOV_SNR"].isna().any():
        raise ValueError("Missing or undefined FOV SNR values")
    print(fovpos.sort_values('FOV_SNR'))

    # Retrieve FOV dimensions (saved by Ashlar or by assembleFOVs.py)
    with open(args.tile_size) as f:
        tile_width, tile_height = map(int, f.read().strip().split(","))
    
    if tile_width <= 0 or tile_height <= 0:
        raise ValueError("FOV dimensions must be positive")

    # Use above SNR values to filter transcripts in overlapping FOV areas
    # (this avoids duplication of transcripts detected in 2 FOVs)
    print("Filtering transcripts in overlapping FOV regions based on pseudo-SNR")
    tx = filter_overlap_tx(tx, fovpos, tile_width=tile_width, tile_height=tile_height)

    # Define indices of transcripts for main, control, and other categories based on regex patterns
    print("Splitting transcripts file by regex patterns for main, control, and other probes")
    idx_ctrl = (
        tx["target"].str.contains(args.ctrl_regex, na=False)
        if args.ctrl_regex
        else pd.Series(False, index=tx.index)
    )
    idx_other = (
        tx["target"].str.contains(args.other_regex, na=False)
        if args.other_regex
        else pd.Series(False, index=tx.index)
    )
    idx_main = ~idx_ctrl & ~idx_other

    # Create PointsModel objects for main and control transcripts
    points = {
        "main": spd.models.PointsModel.parse(
            tx.loc[idx_main], feature_key="target",
            coordinates={"x": "x", "y": "y"}
        ),
        "ctrl": spd.models.PointsModel.parse(
            tx.loc[idx_ctrl], feature_key="target",
            coordinates={"x": "x", "y": "y"}
        ),
    }

    # Only create the 'other' element when a regex is specified
    if args.other_regex:
        points["other"] = spd.models.PointsModel.parse(
            tx.loc[idx_other], feature_key="target",
            coordinates={"x": "x", "y": "y"}
        )

    # If a polygon CSV is provided, read it and create a ShapesModel object
    if args.pg_csv is not None:
        print("Importing polygon CSV file")
        pg = pd.read_csv(args.pg_csv, header=0)
        if "fov" in pg.columns:
            pg.rename(columns={"fov": "FOV"}, inplace=True)
        missing = {"FOV", "cell", "x_global_px", "y_global_px"} - set(pg.columns)
        if missing:
            raise ValueError(f"missing {', '.join(sorted(missing))} in pg table")
        pg = pg.loc[pg["FOV"].isin(fov2pos.keys())].copy()
        pg = transform_coordinates(pg, fov2pos)
        segments = []
        for cell, sub in pg.groupby("cell"):
            coords = list(zip(sub["x"], sub["y"]))
            if len(coords) > 2:
                poly = Polygon(coords)
                segments.append({"cell_id": cell, "geometry": poly})
        nsk = pg["cell"].nunique() - len(segments)
        if nsk > 0:
            print(f"{YELLOW}Warning: {nsk} segments have insufficient points.{RESET}")
        gdf = gpd.GeoDataFrame(segments, geometry="geometry", crs=None)
        gdf.set_index("cell_id", inplace=True)
        shapes = {"atomx_segments": spd.models.ShapesModel.parse(gdf)}
    else:
        shapes = None

    spd.SpatialData(images=images, points=points, shapes=shapes).write(args.output)

    print(f"{GREEN}Zarr dataset created successfully.{RESET}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        die("Failed to create Zarr dataset.")
