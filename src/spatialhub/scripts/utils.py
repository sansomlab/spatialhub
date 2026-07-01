import sys
import traceback

RESET = "\033[0m"
RED = "\033[91m"
GREEN = "\033[92m"
YELLOW = "\033[93m"


def die(msg: str, code: int = 1):
    traceback.print_exc()
    print(f"{RED}{msg}{RESET}")
    sys.exit(code)


def print_arguments(cmdargs):
    print(f"\n{GREEN}Parsed arguments:{RESET}")
    for k, v in vars(cmdargs).items():
        print(f"\t{GREEN}- {k}:\t{v}{RESET}")
    print("")
