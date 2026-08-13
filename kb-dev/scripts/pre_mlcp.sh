#!/bin/bash

set -euo pipefail
source /workspace/system.vars
cd /workspace/simulations/${SYS_NAME}/mlcp

# -- process intermediaries and run MLCP
echo "Running MLCP..."