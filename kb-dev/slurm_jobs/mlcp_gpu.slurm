#!/bin/bash

#SBATCH --job-name=mlcp_gpu
#SBATCH --constraint=gpu         
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --gpus-per-node=4
#SBATCH --time=00:30:00
#SBATCH --output=/global/cfs/cdirs/m5128/kbilal/interface/gpu_job_%j.out
#SBATCH --error=/global/cfs/cdirs/m5128/kbilal/interface/gpu_job_%j.err

# OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

#run the application: 
# For large systems such as naphthalene, the total nr. of GPUs should equal the number of states in top layer node
# The number of GPU nodes to request is (nr. states)/4
echo "Running MLCP (GPU)..."
source system.vars
srun -n 1 -c 128 --gpu-bind=none --cpu_bind=cores podman-hpc run --userns=keep-id --rm --gpu -v "$PWD":/workspace --entrypoint="" working mlcp.x mlcp_${SYS_NAME}.inp > mlcp_${SYS_NAME}.out