#!/bin/bash

set -euo pipefail
source /workspace/system.vars
cd /workspace/simulations/${SYS_NAME}/nwchem

# -- run NWChem
echo "Running NWChem..."
nwchem ${SYS_NAME}.nw > nwc_${SYS_NAME}.out
conda run -n base /workspace/scripts/nwchem_fc.py nwc_${SYS_NAME}.out ${SYS_NAME}

cd /workspace/simulations/${SYS_NAME}/mlcp
cp /workspace/inputs/pbqff2_nmodes_h2o.dat ../pbqff # temporary, until isolated from pbqff
conda run -n base /workspace/scripts/translation.py ${SYS_NAME}
mv ../nwchem/nwcf2${SYS_NAME}.dat f2${SYS_NAME}.dat
mv ../pbqff/f3${SYS_NAME}.dat .
mv ../pbqff/f4${SYS_NAME}.dat .

echo "Translation and NWChem step complete."