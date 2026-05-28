import argparse
import subprocess

from pathlib import Path
from spatialhub import __version__


def welcome_message():
    width = 68

    print()
    print("=" * width)
    print(f"Welcome to SpatialHub (v{__version__})!".center(width))
    print("An Integrated Platform for Spatial Transcriptomics Analysis.".center(width))
    print("=" * width)
    print()


def main():
    welcome_message()

    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="An integrated platform for spatial transcriptomics analysis."
    )
    parser.add_argument("module", help="The module to run.")
    parser.add_argument("task", help="The task to run, [config|full|<rulename>].")
    parser.add_argument("--dry", action="store_true", help="Perform a dry run only.")
    parser.add_argument(
        "--scheduler", default="drmaa", help="The scheduler to use, [drmaa|slurm]."
    )
    parser.add_argument(
        "--cores", default="all", help="Number of cores to use for Snakemake."
    )
    parser.add_argument(
        "--jobs", default="1", help="Number of jobs to run in parallel for Snakemake."
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    args = parser.parse_args()

    # Determine the path to the Snakefile
    snakefile_path = Path(__file__).parent.parent / "workflow" / "Snakefile"

    # Build the Snakemake command
    cmd = [
        "snakemake",
        "--snakefile",
        str(snakefile_path),
        "--config",
        f"module={args.module}",
        f"task={args.task}",
        "--cores",
        str(args.cores),
        "--jobs",
        str(args.jobs),
    ]
    if args.dry:
        cmd.append("--dry-run")
    elif args.scheduler == "drmaa":
        cmd.extend(["--executor", "drmaa"])
        drmaa_args = (
            " -p {resources.partition} --mem={resources.mem_mb} --cpus-per-task={threads} "
            + "--time={resources.time} --output=logs/job_%j.out --error=logs/job_%j.err"
        )
        log_dir = "logs"
        cmd.extend(["--drmaa-args", drmaa_args, "--drmaa-log-dir", log_dir])
    elif args.scheduler == "slurm":
        cmd.extend(["--executor", "slurm"])
    else:
        raise ValueError(f"Unknown scheduler {args.scheduler}, use 'drmaa' or 'slurm'.")

    # Execute the command
    print("Running command:", " ".join(cmd))
    subprocess.run(cmd, check=True)
