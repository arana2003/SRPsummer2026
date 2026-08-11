#!/bin/bash

set -euo pipefail

# Load user-defined system vars
touch system.vars
SYSTEM_VARS_PATH="$(realpath system.vars)"
cat > system.vars << 'EOF'
export WORKING="/global/cfs/cdirs/m5128/kbilal/interface"
export MLCPPIPE="/pscratch/sd/k/kbilal/interface"
export PROJ="m5128"
export EXTEND_PES="false"
export SYS_NAME="methanol"
export POLL_INTERVAL="10"
export INPUT="${MLCPPIPE}/inputs/methanol.inp"
export CHECKPOINT_FILE="${MLCPPIPE}/simulations/${SYS_NAME}/pipeline.checkpoint"
EOF
source system.vars

cleanup() {
    echo "Cleaning up..."
    cp -R ${MLCPPIPE}/simulations/* ${WORKING}/simulations/
    echo "MLCP Pipe closing"
}
trap cleanup EXIT

# Read last completed step (0 if none)
if [ -f "$CHECKPOINT_FILE" ]; then
    LAST_STEP=$(cat "$CHECKPOINT_FILE")
else
    LAST_STEP=0
fi

run_step() {
    STEP_NUM=$1
    STEP_NAME=$2
    shift 2

    if [ "$STEP_NUM" -lt "$LAST_STEP" ]; then
        echo "Skipping step $STEP_NUM ($STEP_NAME) — already completed"
        return 0
    fi

    echo "Running step $STEP_NUM: $STEP_NAME"
    if "$@"; then
        pwd
        echo "$STEP_NUM" > "$CHECKPOINT_FILE"
        echo "Step $STEP_NUM ($STEP_NAME) submitted, checkpoint saved"
    else
        echo "Step $STEP_NUM ($STEP_NAME) failed, stopping"
        return 1
    fi
}

wait_for_job() {
    jobid=$1
    while squeue -h -j "$jobid" 2>/dev/null | grep -q .; do
        sleep "$POLL_INTERVAL"
    done

    sleep 5

    state=$(sacct -j "$jobid" --format=State --noheader --parsable2 | head -1)
    state=${state%% *}      

    if [ "$state" == "FAILED" ]; then
        echo "Job $jobid finished with state: $state" >&2
        return 1
    fi
}

step1_setup() {
    mkdir -p ${MLCPPIPE}/simulations/${SYS_NAME}
    cp ${SYSTEM_VARS_PATH} ${MLCPPIPE}/simulations/${SYS_NAME}/system.vars
    rm -f SYSTEM_VARS_PATH
    SYSTEM_VARS_PATH="${MLCPPIPE}/simulations/${SYS_NAME}/system.vars"
    source ${SYSTEM_VARS_PATH}

    cp -R ${WORKING}/* ${MLCPPIPE}/
    cd ${MLCPPIPE}

    scripts/setup.sh
}
step2_submit_pbqff() {

    echo "Submitting pbqff..."
    source ${SYSTEM_VARS_PATH}
    pbqff_jobid=$(sbatch -A "$PROJ" --export=NONE --parsable slurm_jobs/pbqff.slurm)
    echo "PBQFF_JOBID=$pbqff_jobid" >> "simulations/${SYS_NAME}/jobids"
    echo 'export JOBIDS="simulations/${SYS_NAME}/jobids"' >> ${SYSTEM_VARS_PATH}

    wait_for_job "$pbqff_jobid"
}

step3_submit_nwchem() {

    echo "Submitting nwchem..."
    source ${SYSTEM_VARS_PATH}
    nwchem_jobid=$(command sbatch -A ${PROJ} --parsable --dependency=afterok:$pbqff_jobid --kill-on-invalid-dep=yes slurm_jobs/nwchem.slurm)
    echo "NWCHEM_JOBID=$nwchem_jobid" >> "simulations/${SYS_NAME}/jobids"

    wait_for_job "$nwchem_jobid"
}

step3_submit_extension() {

    echo "Submitting extension..."
    source ${SYSTEM_VARS_PATH}
    extension_jobid=$(command sbatch -A ${PROJ} --parsable --dependency=afterok:$pbqff_jobid --kill-on-invalid-dep=yes slurm_jobs/extension.slurm)
    echo "EXTENSION_JOBID=$extension_jobid" >> "simulations/${SYS_NAME}/jobids"

    wait_for_job "$extension_jobid"
}

step4_submit_mlcp() {
    echo "Submitting mlcp (cpu)..."
    source ${SYSTEM_VARS_PATH}

    # mlcp_cpu_jobid=$(command sbatch -A ${PROJ} --parsable --dependency=afterok:$nwchem_jobid --kill-on-invalid-dep=yes slurm_jobs/mlcp_cpu.slurm)
    # echo "MLCP_CPU_JOBID=$mlcp_cpu_jobid" >> "simulations/${SYS_NAME}/jobids"

    # wait_for_job "$mlcp_cpu_jobid"

    # echo "Submitting mlcp (gpu)..."
    # sed -i 's/donode="F"/donode="T"/g' simulations/${SYS_NAME}/mlcp/mlcp_${SYS_NAME}.inp
    mlcp_gpu_jobid=$(command sbatch -A ${PROJ} --parsable --dependency=afterok:$nwchem_jobid --kill-on-invalid-dep=yes slurm_jobs/mlcp_gpu.slurm)
    echo "MLCP_GPU_JOBID=$mlcp_gpu_jobid" >> "simulations/${SYS_NAME}/jobids"

    wait_for_job "$mlcp_gpu_jobid"
}

if [[ "${1:-}" == "-r" || "${1:-}" == "--reset" ]]; then
    LAST_STEP=0
fi

run_step 1 "setup" step1_setup

run_step 2 "pbqff" step2_submit_pbqff

if [ "$EXTEND_PES" = false ]; then
    podman-hpc run --userns=keep-id --rm --gpu -v "$PWD":/workspace --entrypoint="/workspace/scripts/pre_nwchem.sh" working 
    run_step 3 "nwchem" step3_submit_nwchem
    podman-hpc run --userns=keep-id --rm --gpu -v "$PWD":/workspace --entrypoint="/workspace/scripts/post_nwchem.sh" working 
else
    run_step 3 "extension" step3_submit_extension
fi

podman-hpc run --userns=keep-id --rm --gpu -v "$PWD":/workspace --entrypoint="/workspace/scripts/pre_mlcp.sh" working 
run_step 4 "mlcp" step4_submit_mlcp

echo "GPU step complete."
echo 'Completed MLCP!'

echo "Pipeline complete!"