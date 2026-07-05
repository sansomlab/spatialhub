import subprocess

from argparse import ArgumentParser as AP
from importlib.resources import files
from pathlib import Path
from spatialhub import __version__

RED = "\033[91m"
GREEN = "\033[92m"
RESET = "\033[0m"


def welcome_message():
    width = 68

    print()
    print("=" * width)
    print(f"Welcome to SpatialHub (v{__version__})!".center(width))
    print("Pipelines for Spatial Transcriptomics Analysis.".center(width))
    print("=" * width)
    print()


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


def run_snakemake(smkpath, task, cores, jobs, dry_run, lock):
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
    cmd.extend(["--executor", "drmaa"])
    drmaa_args = (
        " -A {resources.account} -p {resources.partition} "
        + "--mem={resources.mem_mb} --cpus-per-task={threads} --time={resources.time} "
        + "--output=logs/job_%j.out --error=logs/job_%j.err "
    )
    log_dir = "logs"
    cmd.extend(["--drmaa-args", drmaa_args, "--drmaa-log-dir", log_dir])
    if dry_run:
        cmd.append("--dry-run")
    cmd.extend([task])
    print("Running command:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    welcome_message()

    vinfo = f"%(prog)s {__version__}"

    # Parse command-line arguments
    p = AP(description="Pipelines for Spatial Transcriptomics Analysis.")
    p.add_argument("workflow", help="The workflow to run.")
    p.add_argument("task", help="The task to run, [config|full|<rulename>].")
    p.add_argument("--dry", action="store_true", help="Perform a dry run only.")
    p.add_argument("--cores", default="all", help="Number of cores to use.")
    p.add_argument("--jobs", default="1", help="Number of parallel jobs to run.")
    p.add_argument("--lock", action="store_true", help="Lock output by chmod.")
    p.add_argument("-v", "--version", action="version", version=vinfo)
    args = p.parse_args()

    # Determine the path to the Snakefile
    smkpath = files("spatialhub").joinpath("workflow", f"{args.workflow}.smk")

    if args.lock and args.task == "config":
        print(f"{GREEN}[INFO] `--lock` is ignored when running config task.{RESET}")

    if args.task == "config":
        copy_config_file(args.workflow)
    else:
        run_snakemake(
            smkpath,
            args.task,
            str(args.cores),
            str(args.jobs),
            args.dry,
            str(args.lock).lower(),
        )

    print(f"{GREEN}Completed successfully!{RESET}")
