#!/bin/bash
#SBATCH --job-name=nwchem
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

# On one Perlmutter GPU node, make sure the arguments of -n and -c multiply to 128 with -c larger than 1 (to avoid running on hyperthreads), so you can use anything from -n 1 -c 128 to -n 64 -c 2 .

# Old call: srun -n 1 -c 128 --cpu_bind=cores -G 4 --gpu-bind=none podman-hpc run --userns=keep-id --rm --gpu -v "$PWD":/workspace --entrypoint="/workspace/scripts/nwchem.sh" working
source system.vars
srun -n 1 -c 128 --gpu-bind=none --cpu_bind=cores podman-hpc run --userns=keep-id --rm --gpu -v "$PWD":/workspace --entrypoint="" working nwchem ${SYS_NAME}.nw > nwc_${SYS_NAME}.out