#!/usr/bin/env bash
set -euo pipefail

# local install for interactive use

micromamba create -f environment.yaml -y
micromamba run -n ngs-analysis python -m ipykernel install --user --name=ngs-analysis
