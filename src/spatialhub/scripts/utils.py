import sys
import traceback

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
