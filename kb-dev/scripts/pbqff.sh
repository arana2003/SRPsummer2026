#!/bin/bash

set -euo pipefail
source /workspace/system.vars
cd /workspace/simulations/${SYS_NAME}/pbqff

# -- run PBQFF
echo "Running PBQFF..."

if [ -f "pbqff.out" ]; then
    pbqff -c ${SYS_NAME}.toml
else
    pbqff -o ${SYS_NAME}.toml
fi

conda run -n base /workspace/scripts/qfflist2.py pbqff.out ${SYS_NAME}

echo "CPU step complete."