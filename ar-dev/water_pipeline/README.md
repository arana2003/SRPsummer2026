# Alexis Rana — Multifidelity PES Pipeline
## SRP 2026 | LBNL AMCRD | PI: Tucker Carrington | Mentor: Dr. Phillip Thomas

Automated multifidelity PES fitting combining MOPAC PM7 (low-fidelity),
NWChem B3LYP/6-31G* (high-fidelity), and PBQFF 3-mode force constants.
Output: f2/f3/f4_mt.dat files for MLCP.

## Quick Start (run from repo root)

### Step 1 - Save PBQFF data
Extract geometry and NC modes from pbqff.out into .npy files (see full README).

### Step 2 - Calibrate scale factors
podman-hpc run --rm -v "$(pwd)":/pipeline -v "$SCRATCH":/scratch -v /dev/shm:/dev/shm working bash -c "
/opt/miniconda3/bin/python /pipeline/ar-dev/water_pipeline/calibrate_modes.py
    --mol water --geom /scratch/water_modes/water_pbqff_geom.npy
    --modes /scratch/water_modes/water_pbqff_NC_modes.npy
    --freqs /scratch/water_modes/water_pbqff_NC_freqs.npy
    --labels /scratch/water_modes/water_pbqff_labels.npy
    --equil-hof -57.79981978 --mopac /opt/mopac/mopac
    --outfile /scratch/water_modes/water_pbqff_scales.npy --shm /dev/shm"

### Step 3 - 1D MOPAC grid (55 pts per mode)
for mode in 1 2 3; do
    podman-hpc run --rm -v "$(pwd)":/pipeline -v "$SCRATCH":/scratch -v /dev/shm:/dev/shm
        -e MOLDATA_DIR=/scratch/water_modes working bash -c "
    /opt/miniconda3/bin/python /pipeline/ar-dev/water_pipeline/worker_v3.py
        --mol water --mode-i $mode --mode-j $mode
        --pointsfile /scratch/points_1d.txt
        --outcsv /scratch/results_1d_mode${mode}.csv
        --equil-hof -57.79981978 --mopac /opt/mopac/mopac --shm /dev/shm" &
done && wait

### Step 4 - 2D MOPAC grid (50 Sobol pts per pair)
python3 ar-dev/water_pipeline/generate_2d_points.py --mol water --n-train 50 --q-max 1.5 --outdir $SCRATCH/mopac_2d_water
podman-hpc run --rm -v "$SCRATCH/mopac_2d_water":/mopac working bash -c "cd /mopac && for f in *.mop; do /opt/mopac/mopac \$f; done"

### Step 5 - Convert PBQFF f3/f4 to Morse/tanh SOP
python3 ar-dev/water_pipeline/f34_to_morse_tanh_sop.py
    --f3 kb-dev/tests/h2o/mlcp/pes/f3h2o.dat
    --f4 kb-dev/tests/h2o/mlcp/pes/f4h2o.dat
    --morse ar-dev/water_pipeline/Morse_water.dat
    --outfile water_sop.dat --system h2o

### Step 6 - Run MLCP (toggle system name to switch PES)
# system='h2o'    -> PBQFF only
# system='h2o_mt' -> Alexis multifidelity pipeline
cp ar-dev/water_pipeline/f2h2o_mt.dat kb-dev/tests/h2o/mlcp/pes/
cp ar-dev/water_pipeline/f3h2o_mt.dat kb-dev/tests/h2o/mlcp/pes/
cp ar-dev/water_pipeline/f4h2o_mt.dat kb-dev/tests/h2o/mlcp/pes/
/path/to/mlcp mlcp_h2o_mt.inp > mlcp_h2o_mt.out

## Key Files
calibrate_modes.py       - Empirical per-mode scale calibration
worker_v3.py             - Persistent MOPAC workers, /dev/shm I/O
generate_2d_points.py    - Sobol 2D grid generation
fit_sklearn_gpr.py       - 1D GPR fitting
fit_2d_multifidelity.py  - 2D delta-GPR fitting
f34_to_morse_tanh_sop.py - PBQFF f3/f4 Taylor to Morse/tanh conversion
mfgp.py                  - Multifidelity GPR implementation
gpr_to_sop_mlcp.py       - GPR to MLCP SOP format
fit_tanh_sop.py          - Tanh SOP fitting
fit_morse_sop.py         - Morse SOP fitting
generate_job_list.py     - Flat job list for worker dispatch
run_parallel.sh          - Parallel run script for HPC
load_data.py             - Data loading utilities
load_test_data.py        - Test set loading utilities
read_kaiwan_outputs.py   - Interface with Kaiwan PBQFF/NWChem pipeline
get_normal_modes.py      - MOPAC OPT+FORCE (cross-check only)
detect_symmetry.py       - Auto-detect symmetric vs asymmetric modes
parse_sample_inp.py      - Parse Kaiwan sample.inp format
regenerate_test_set.py   - Regenerate test set after fixes
fit_multifidelity_combined.py - Combined multifidelity fitting
Morse_water.dat          - Morse parameters (freq, alpha, shift)
f2h2o_mt.dat             - MLCP: 1D multifidelity (tanh basis)
f3h2o_mt.dat             - MLCP: 2D GPR + 3D PBQFF cubic
f4h2o_mt.dat             - MLCP: 2D GPR + 3D PBQFF quartic

## Validated Results (Water)
V(q=0) = 0.00 cm-1 exactly for all 3 modes
V(q=0.3)/V_harm = 0.983, 0.984, 0.967 (real anharmonicity)
2D GPR RMSE: 8.81, 3.79, 5.38 cm-1 for pairs (1,2),(1,3),(2,3)
3-mode coupling y1*y2*y3 = 10.47 cm-1 (small for water)

## Key Technical Notes
- PBQFF geometry in pbqff.out is in Angstrom (NOT Bohr)
- PBQFF NC matrix: columns = mode vectors, highest to lowest freq
- Do NOT use pbqff2_nmodes_{nm}.dat (NWChem vectors, not PBQFF)
- Do NOT use MOPAC FORCE modes (violate molecular symmetry)
- Displacement: delta_x_ang = q * scale * L where scale = sqrt(nu/k)
- No single formula for q conversion - use calibrate_modes.py

## Contact
Alexis Rana | arana2003@gmail.com
OHSU PhD Fall 2026 | SRP 2026 LBNL AMCRD
