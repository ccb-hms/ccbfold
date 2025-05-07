[![License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](./LICENSE)
# CCB/RC Alphafold 3 workshop

The goal of this workshop is to provide an introduction to running Alphafold 3 jobs on the O2 cluster at HMS

These scripts have been tested to be capeable of running on the cluster.


## Prerequisites
- Knowledge of:
    - Running jobs on O2
    - Bash/terminal basics
    - Python virtural environments and package management
- Software installed on your laptop
    - Python >=3.12
    - [uv](https://docs.astral.sh/uv/)
    - [juv](https://github.com/manzt/juv)

## Install

```bash
uv sync
#for exercise05 visualization
uv sync --extra visualization
uv pip install -e .
source .venv/bin/activate
```

# OR if omitting uv
```bash
python3 pip install -e .
```