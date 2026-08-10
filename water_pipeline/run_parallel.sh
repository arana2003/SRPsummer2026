#!/bin/bash
# Run N_WORKERS inside one container using Python multiprocessing
# Usage: bash run_parallel.sh water 1 2 128 1024 -57.01295260828556

MOLECULE=${1:-water}
MODE_I=${2:-1}
MODE_J=${3:-2}
N_WORKERS=${4:-128}
JOBS_PER_NODE=${5:-32000}
EQUIL_HOF=${6:-0.0}
POINTSFILE=${7:-$SCRATCH/points.txt}
OUTCSV=${8:-$SCRATCH/results.csv}

podman-hpc run --rm \
    -v "$HOME/SRPsummer2026":/pipeline \
    -v "$SCRATCH":/scratch \
    -v /dev/shm:/dev/shm \
    working bash -c "
/opt/miniconda3/bin/python - << PYEOF
import subprocess, multiprocessing

def run_worker(worker_id):
    subprocess.run([
        '/opt/miniconda3/bin/python',
        '/pipeline/water_pipeline/worker_v3.py',
        '--worker-id', str(worker_id),
        '--n-workers', '$N_WORKERS',
        '--node-id', '0',
        '--jobs-per-node', '$JOBS_PER_NODE',
        '--mol', '$MOLECULE',
        '--mode-i', '$MODE_I',
        '--mode-j', '$MODE_J',
        '--pointsfile', '/scratch/$(basename $POINTSFILE)',
        '--outcsv', '/scratch/$(basename $OUTCSV)',
        '--equil-hof', '$EQUIL_HOF',
        '--mopac', '/opt/mopac/mopac',
        '--shm', '/dev/shm',
    ])

with multiprocessing.Pool($N_WORKERS) as pool:
    pool.map(run_worker, range($N_WORKERS))
print('Done')
PYEOF
"
echo "Results: \$(wc -l < $OUTCSV) rows"
