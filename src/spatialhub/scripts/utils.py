import sys
import traceback

from spatialdata import bounding_box_query

RESET = "\033[0m"
RED = "\033[91m"
GREEN = "\033[92m"
YELLOW = "\033[93m"


def die(msg: str, code: int = 1):
    assert code != 0, "Exit code must be non-zero for die()"
    traceback.print_exc()
    print(f"{RED}{msg}{RESET}")
    sys.exit(code)


def print_arguments(cmdargs):
    print(f"\n{GREEN}Parsed arguments:{RESET}")
    for k, v in vars(cmdargs).items():
        print(f"\t{GREEN}- {k}:\t{v}{RESET}")
    print("")


def parse_fov_list(fov_lst: str):
    """
    Parse a comma-separated list of FOVs and ranges.

    Args
    ------
        fov_lst: Comma-separated list of FOVs and ranges (e.g. "1,2,5-7,10").

    Returns
    -------
        list[int]: Parsed FOV identifiers.

    Raises
    -------
        ValueError: If the input contains malformed ranges or invalid FOV identifiers.
    """
    fovs = []
    for fov_str in fov_lst.split(","):
        fov_str = fov_str.strip()
        if not fov_str:
            continue

        if "-" in fov_str:
            try:
                start, end = map(int, fov_str.split("-"))
            except ValueError:
                raise ValueError(f"invalid FOV range '{fov_str}'")
            if start > end:
                raise ValueError(f"invalid FOV range '{fov_str}': start > end")
            fovs.extend(range(start, end + 1))
        else:
            try:
                fovs.append(int(fov_str))
            except ValueError:
                raise ValueError(f"invalid FOV '{fov_str}'")
    return fovs


def crop_sdata(sdata, axes, x_minmax, y_minmax, coords="global"):
    """
    Crop a SpatialData object to a specified bounding box.

    Args
    ------
        sdata: SpatialData object to crop.
        axes: List of axes to consider.
        x_minmax: List containing the minimum and maximum x coordinates of the bounding box.
        y_minmax: List containing the minimum and maximum y coordinates of the bounding box.
        coords: Coordinate system for the bounding box.

    Returns
    -------
        SpatialData: Cropped SpatialData object.
    """
    return bounding_box_query(
        sdata,
        axes=axes,
        min_coordinate=[x_minmax[0], y_minmax[0]],
        max_coordinate=[x_minmax[1], y_minmax[1]],
        target_coordinate_system=coords,
    )
