#!/bin/bash

set -euo pipefail
source /workspace/system.vars
cd /workspace/simulations/${SYS_NAME}/mlcp

# -- process intermediaries and run MLCP
echo "Running MLCP..."
mlcp.x mlcp_${SYS_NAME}.inp > mlcp_${SYS_NAME}.out

echo "GPU step complete."
echo 'Completed MLCP pipeline!'