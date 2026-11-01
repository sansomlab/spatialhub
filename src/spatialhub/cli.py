import subprocess
import shutil

import yaml

from argparse import ArgumentParser as AP
from importlib.resources import files
from pathlib import Path

from rich_argparse import RawDescriptionRichHelpFormatter

from spatialhub import __version__

RED = "\033[91m"
GREEN = "\033[92m"
YELLOW = "\033[93m"
RESET = "\033[0m"

PLACEHOLDER_ACCOUNT = "slurm_account"

# Executors whose configs use native Slurm resource names.
NATIVE_NAME_EXECUTORS = ("slurm-drmaa", "slurm")

# Single source of truth for available workflows.
#
# executor:
#   "drmaa" - submitted through the DRMAA executor (CPU workflows).
#   "slurm-drmaa" - submitted through snakemake-executor-plugin-slurm-drmaa,
#             which builds job options from each rule's resources and supports
#             GPU requests.
#   "slurm" - submitted through the native Slurm executor plugin.
#   Configs for "slurm-drmaa" and "slurm" workflows use native resource names:
#   slurm_account, slurm_partition, runtime (minutes) and gres.
WORKFLOWS = {
    "cosmx_makeImage": {
        "summary": "Generate one image file per CosMx sample from raw FOV data.",
        "steps": [
            "with Ashlar (stitching FOVs): genBlankFOV -> completeGrid -> runAshlar",
            "without Ashlar (appending FOVs without resolving overlaps at "
            "boundaries): assembleFOVs",
        ],
        "executor": "drmaa",
    },
    "cosmx_makeZarr": {
        "summary": "Generate one Zarr file per CosMx sample from compiled image "
        "and raw transcripts data.",
        "steps": [],
        "executor": "drmaa",
    },
    "visiumhd_makeZarr": {
        "summary": "Generate one Zarr file per Visium HD capture area from FASTQ "
        "and images.",
        "steps": ["Space Ranger count", "Zarr generation"],
        "executor": "drmaa",
    },
    "visiumhd_prepareBin2Cell": {
        "summary": "Generate intermediate files up to the StarDist step of Bin2Cell.",
        "steps": [],
        "executor": "drmaa",
    },
    "extractH5AD": {
        "summary": "Build an AnnData (H5AD) from a SpatialData (Zarr) when a table "
        "already exists (typically Visium HD).",
        "steps": [],
        "executor": "drmaa",
    },
    "aggregateH5AD": {
        "summary": "Build an AnnData (H5AD) from a SpatialData (Zarr) when no table "
        "exists, by aggregating points and shapes (typically CosMx).",
        "steps": [],
        "executor": "drmaa",
    },
    "addShapes": {
        "summary": "Add shape layers (Cellpose cell segmentation, FOV boxes) to "
        "existing SpatialData (Zarr) files. Uses GPUs.",
        "steps": [],
        "executor": "slurm-drmaa",
    },
    "runRCTD": {
        "summary": "Deconvolve cell-type abundances against a single-cell reference "
        "using RCTD.",
        "steps": [],
        "executor": "drmaa",
    },
    "runCell2Location": {
        "summary": "Deconvolve cell-type abundances against a single-cell reference "
        "using Cell2Location.",
        "steps": [],
        "executor": "drmaa",
    },
}


def workflows_help():
    """Render a compact workflow catalogue for the --help epilog.

    Uses rich markup (interpreted by RawDescriptionRichHelpFormatter). Only the
    one-line summary is shown per workflow; per-step detail is intentionally
    omitted here to keep the top-level help scannable.
    """
    width = max(len(name) for name in WORKFLOWS)
    lines = ["[bold]Available workflows:[/]", ""]
    for name, info in WORKFLOWS.items():
        lines.append(f"  [bold cyan]{name:<{width}}[/]  {info['summary']}")
    lines.append("")
    lines.append("Run a workflow with:  [bold]spatialhub <workflow> <task>[/]")
    lines.append("  <task> is one of: config | full | <rulename>")
    return "\n".join(lines)


def welcome_message():
    width = 68

    print()
    print("=" * width)
    print(f"Welcome to SpatialHub (v{__version__})!".center(width))
    print("Workflows for Spatial Transcriptomics Analysis.".center(width))
    print("=" * width)
    print()


def get_default_slurm_account():
    """Return the user's default Slurm account, or None if it cannot be determined."""
    if shutil.which("sacctmgr") is None:
        return None
    try:
        result = subprocess.run(
            ["sacctmgr", "-nP", "show", "user", "$USER", "format=DefaultAccount"],
            capture_output=True,
            text=True,
            timeout=15,
        )
    except (subprocess.SubprocessError, OSError):
        return None
    if result.returncode != 0:
        return None
    account = result.stdout.strip()
    return account or None


def get_valid_slurm_accounts():
    """Return the set of accounts the user is associated with, or None if undeterminable."""
    if shutil.which("sacctmgr") is None:
        return None
    try:
        result = subprocess.run(
            ["sacctmgr", "-nP", "show", "assoc", "user=$USER", "format=Account"],
            capture_output=True,
            text=True,
            timeout=15,
        )
    except (subprocess.SubprocessError, OSError):
        return None
    if result.returncode != 0:
        return None
    accounts = {line.strip() for line in result.stdout.splitlines() if line.strip()}
    return accounts or None


def _config_error(message):
    print(f"{RED}[ERROR] {message}{RESET}")
    raise SystemExit(1)


def _account_error(message):
    print(f"{RED}[ERROR] {message}{RESET}")
    default_account = get_default_slurm_account()
    if default_account:
        print(
            f"{YELLOW}[HINT] Your default Slurm account appears to be "
            f"'{default_account}'. Set 'resources.slurm_account' in the workflow "
            f"config accordingly.{RESET}"
        )
    else:
        print(
            f"{YELLOW}[HINT] To find your default Slurm account, run:\n"
            f"           sacctmgr -nP show user $USER format=DefaultAccount{RESET}"
        )
    raise SystemExit(1)


def validate_config(workflow):
    """Pre-flight check of the workflow config before Snakemake is invoked.

    Surfaces common configuration mistakes here rather than inside a verbose
    Snakemake traceback:
      - the Slurm account is missing or still the template placeholder;
      - for native-Slurm workflows, legacy key names that the native executor
        silently ignores ('account' instead of 'slurm_account', 'time' instead
        of 'runtime').
    """
    executor = WORKFLOWS[workflow]["executor"]
    config_path = Path.cwd().joinpath(f"{workflow}.yaml")
    if not config_path.exists():
        _config_error(
            f"Config file '{config_path.name}' not found in the current directory. "
            f"Run 'spatialhub {workflow} config' to create it."
        )

    with config_path.open() as fh:
        config = yaml.safe_load(fh) or {}

    resources = config.get("resources") or {}

    if executor in NATIVE_NAME_EXECUTORS:
        if not resources.get("slurm_account") and resources.get("account"):
            _config_error(
                f"'{workflow}' runs on the {executor} executor, which reads "
                f"'slurm_account', not 'account'. Rename the key under 'resources' "
                f"in {config_path.name}."
            )
        account = resources.get("slurm_account")
    else:
        account = resources.get("slurm_account") or resources.get("account")

    if not account:
        _account_error(
            f"No Slurm account set under 'resources.slurm_account' in "
            f"{config_path.name}."
        )

    if account == PLACEHOLDER_ACCOUNT:
        _account_error(
            f"The Slurm account in {config_path.name} is still the template "
            f"placeholder '{PLACEHOLDER_ACCOUNT}'. Assign a valid Slurm account."
        )

    if executor in NATIVE_NAME_EXECUTORS:
        for section in ("resources", "gpu_resources"):
            block = config.get(section) or {}
            if "time" in block and "runtime" not in block:
                _config_error(
                    f"'{section}.time' in {config_path.name} is not read by the "
                    f"{executor} executor. Replace it with 'runtime' in minutes "
                    f"(for example, runtime: 240 for 4 hours)."
                )

    valid_accounts = get_valid_slurm_accounts()
    if valid_accounts is not None and account not in valid_accounts:
        _account_error(
            f"Slurm account '{account}' in {config_path.name} is not one of "
            f"your associated accounts ({', '.join(sorted(valid_accounts))})."
        )

    print(f"{GREEN}[INFO] Config validated: using Slurm account '{account}'.{RESET}")


def copy_config_file(module):
    yaml_dir = Path(__file__).parent.joinpath("configfiles")
    src = yaml_dir.joinpath(f"{module}.yaml")
    dst = Path.cwd().joinpath(f"{module}.yaml")
    if not dst.exists():
        dst.write_text(src.read_text())
        print(f"[INFO] Copied {module} config file to the current directory.")
    else:
        print(f"{RED}[ERROR] {module} configfile already exists.{RESET}")
        raise FileExistsError()


def run_snakemake(smkpath, task, cores, jobs, dry_run, lock, executor):
    cmd = [
        "snakemake",
        "--snakefile",
        str(smkpath),
        "--cores",
        cores,
        "--jobs",
        jobs,
        "-p",
    ]
    cmd.extend(["--config", f"lock={str(lock).lower()}"])

    if executor == "drmaa":
        cmd.extend(["--executor", "drmaa"])
        drmaa_args = (
            " -A {resources.account} -p {resources.partition} "
            + "--mem={resources.mem_mb} --cpus-per-task={threads} --time={resources.time} "
            + "--output=logs/job_%j.out --error=logs/job_%j.err "
        )
        log_dir = "logs"
        cmd.extend(["--drmaa-args", drmaa_args, "--drmaa-log-dir", log_dir])
    elif executor == "slurm-drmaa":
        # The slurm-drmaa executor builds each job's options from its resources
        # (slurm_account, slurm_partition, runtime, mem_mb, threads, gres).
        cmd.extend(["--executor", "slurm-drmaa", "--slurm-drmaa-log-dir", "logs"])
    elif executor == "slurm":
        # The native executor reads slurm_account, slurm_partition, runtime,
        # mem_mb, threads and gres from each rule's resources directly, and
        # writes job logs under .snakemake/slurm_logs/.
        cmd.extend(["--executor", "slurm"])
    else:
        raise ValueError(f"Unknown executor '{executor}'.")

    if dry_run:
        cmd.append("--dry-run")
    cmd.extend([task])
    print("Running command:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    welcome_message()

    vinfo = f"%(prog)s {__version__}"

    # Parse command-line arguments
    p = AP(
        description="Workflows for Spatial Transcriptomics Analysis.",
        epilog=workflows_help(),
        formatter_class=RawDescriptionRichHelpFormatter,
    )
    p.add_argument(
        "workflow",
        choices=WORKFLOWS.keys(),
        metavar="workflow",
        help="The workflow to run (see the list below).",
    )
    p.add_argument("task", help="The task to run, [config|full|<rulename>].")
    p.add_argument("--dry", action="store_true", help="Perform a dry run only.")
    p.add_argument("--cores", default="all", help="Number of cores to use.")
    p.add_argument("--jobs", default="1", help="Number of parallel jobs to run.")
    p.add_argument("--lock", action="store_true", help="Lock output by chmod.")
    p.add_argument("-v", "--version", action="version", version=vinfo)
    args = p.parse_args()

    # Determine the path to the Snakefile
    smkpath = files("spatialhub").joinpath("workflow", f"{args.workflow}.smk")
    executor = WORKFLOWS[args.workflow]["executor"]

    if args.lock and args.task == "config":
        print(f"{GREEN}[INFO] `--lock` is ignored when running config task.{RESET}")

    if args.task == "config":
        copy_config_file(args.workflow)
    else:
        validate_config(args.workflow)
        run_snakemake(
            smkpath,
            args.task,
            str(args.cores),
            str(args.jobs),
            args.dry,
            str(args.lock).lower(),
            executor,
        )

    print(f"{GREEN}Completed successfully!{RESET}")
