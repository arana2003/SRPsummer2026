#!/bin/bash

set -euo pipefail
source /workspace/inputs/system.vars
cd /workspace/simulations/${SYS_NAME}/nwchem

# -- run NWChem
echo "Running NWChem..."