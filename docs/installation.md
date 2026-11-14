# Installation

<div class="nord" markdown="1">
```py
uv venv -p 3.12 spatialhub-py312
source spatialhub-py312/bin/activate
python --version  # Python 3.12.12

# Clone the repo and change directory to it 
git clone https://github.com/sansomlab/spatialhub.git
cd spatialhub     # OR /path/to/clone/repo 

uv pip install .
spatialhub --version
```


## Usage

SpatialHub provides the following command-line interface (CLI). A minimal example is shown below:


```py
spatialhub $workflow_name full
```

For the available `$workflow_name` values, please refer to the [Workflow⁠](#workflows) section below.

The complete command-line interface is shown below:

```py
usage: spatialhub [-h] [--dry] [--cores CORES] [--jobs JOBS] [--lock] [-v] workflow task

Pipelines for Spatial Transcriptomics Analysis.

positional arguments:
  workflow       The workflow to run.
  task           The task to run, [config|full|<rulename>].

options:
  -h, --help     show this help message and exit
  --dry          Perform a dry run only.
  --cores CORES  Number of cores to use.
  --jobs JOBS    Number of parallel jobs to run.
  --lock         Lock output by chmod.
  -v, --version  show program's version number and exit
```

A special feature is provided via the `--lock` option. When enabled, SpatialHub applies `chmod -R a-w` to the output directory to prevent accidental modification or deletion of the results.

To unlock the output directory for modification or deletion, run:

```py
chmod -R +w /folder/to/unlock
```

A built-in task, `config`, is provided to generate a template configuration YAML file, i.e.,

```py
spatialhub $workflow_name config
```
</div>