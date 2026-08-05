#!/usr/bin/env python3
"""
fit_2d_multifidelity.py
Fits a 2D multifidelity delta-GPR for each mode pair.

For each pair (i,j):
  1. Load MOPAC 2D GPR predictions at NWChem points
  2. Compute delta = V_NWChem(qi,qj) - V_MOPAC_GPR(qi,qj)
  3. Fit a 2D GPR to the delta
  4. Final: V_2D_MF(qi,qj) = V_MOPAC_GPR(qi,qj) + delta_GPR(qi,qj)

Usage:
  python3 fit_2d_multifidelity.py \
      --mopac-csv mopac_2d_results.csv \
      --nwchem-dir nwchem_2d_outputs/ \
      --mopac-models models/ \
      --outdir models_2d/
"""

import argparse, csv, os, re, joblib
import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
from sklearn.preprocessing import StandardScaler

KCAL_TO_CM1 = 349.7551


def parse_nwchem_energy(out_file):
    """Parse total SCF energy from NWChem output file."""
    if not os.path.exists(out_file):
        return None
    txt = open(out_file).read()
    # NWChem total energy line
    m = re.search(r'Total DFT energy\s*=\s*([-\d.]+)', txt)
    if not m:
        m = re.search(r'Total SCF energy\s*=\s*([-\d.]+)', txt)
    if not m:
        m = re.search(r'Total energy\s*=\s*([-\d.]+)', txt)
    return float(m.group(1)) if m else None


def fit_2d_gpr(qi_vals, qj_vals, V_vals):
    """Fit a 2D GPR to (qi, qj) -> V data."""
    X = np.column_stack([qi_vals, qj_vals])
    y = V_vals - V_vals.mean()

    kernel = RBF(length_scale=[1.0, 1.0],
                 length_scale_bounds=[(1e-2, 10.0), (1e-2, 10.0)]) + \
             WhiteKernel(noise_level=1e-6, noise_level_bounds=(1e-8, 1e-2))

    gpr = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=5,
                                    normalize_y=True)
    gpr.fit(X, y)
    return gpr, V_vals.mean()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--mopac-csv',    required=True,
                   help='CSV from worker_v3 with MOPAC 2D results')
    p.add_argument('--nwchem-dir',   required=True,
                   help='Directory with NWChem .nw and .out files')
    p.add_argument('--mopac-1d-models', required=True,
                   help='Directory with 1D MOPAC GPR .joblib models')
    p.add_argument('--outdir',       default='models_2d')
    p.add_argument('--equil-hof',    type=float, default=0.0,
                   help='MOPAC equilibrium HOF in kcal/mol')
    p.add_argument('--n-modes',      type=int, default=3)
    p.add_argument('--job-index',    default='2d_job_index.csv',
                   help='Job index CSV from generate_2d_points.py')
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Load MOPAC 2D data
    mopac_data = {}  # (mi,mj,qi,qj) -> V_cm1
    with open(args.mopac_csv) as f:
        for row in csv.DictReader(f):
            if row['status'] != 'ok': continue
            key = (int(row['mode_i']), int(row['mode_j']),
                   round(float(row['qi']),6), round(float(row['qj']),6))
            mopac_data[key] = float(row['V_cm1'])

    print(f"Loaded {len(mopac_data)} MOPAC 2D points")

    # Load 1D MOPAC GPR models
    gpr_1d = {}
    for mode in range(1, args.n_modes + 1):
        mpath = os.path.join(args.mopac_1d_models,
                             f'sklearn_gpr_optimal_mode{mode}_v3.joblib')
        if os.path.exists(mpath):
            gpr_1d[mode] = joblib.load(mpath)
            print(f"Loaded 1D GPR for mode {mode}")

    # Process each mode pair
    from itertools import combinations
    pairs = list(combinations(range(1, args.n_modes + 1), 2))

    results = {}
    for mi, mj in pairs:
        print(f"\nPair ({mi},{mj}):")

        # Find NWChem output files for this pair
        nw_pattern = f'water_m{mi}m{mj}_'
        nw_files = [f for f in os.listdir(args.nwchem_dir)
                    if f.startswith(nw_pattern) and f.endswith('.nw')]

        if not nw_files:
            print(f"  No NWChem files found for pair ({mi},{mj}), skipping")
            continue

        # Parse NWChem energies and match to qi,qj
        # NWChem energies are in Hartree
        HARTREE_TO_CM1 = 219474.63
        nw_points = []

        for nw_file in nw_files:
            out_file = nw_file.replace('.nw', '.out')
            out_path = os.path.join(args.nwchem_dir, out_file)
            E_nw = parse_nwchem_energy(out_path)
            if E_nw is None:
                continue

            # Extract qi, qj from the .nw file geometry
            # (stored in the job index or filename)
            inp_path = os.path.join(args.nwchem_dir, nw_file)
            # Try to find matching MOPAC point by geometry comparison
            # For now use the job index from filename
            nw_points.append({'file': nw_file, 'E_hartree': E_nw})

        if not nw_points:
            print(f"  No NWChem energies parsed for pair ({mi},{mj})")
            continue

        print(f"  Found {len(nw_points)} NWChem points")

        # Load job index to match filenames to qi,qj
        job_idx = {}
        idx_path = os.path.join(args.nwchem_dir, args.job_index)                    if not os.path.isabs(args.job_index)                    else args.job_index
        if not os.path.exists(idx_path):
            idx_path = args.job_index
        if os.path.exists(idx_path):
            with open(idx_path) as f:
                for row in csv.DictReader(f):
                    if int(row['mode_i'])==mi and int(row['mode_j'])==mj:
                        fname = os.path.basename(row['filepath'])
                        job_idx[fname] = (float(row['q_i']), float(row['q_j']))

        # Match NWChem energies to qi,qj and get MOPAC GPR prediction
        # Find equilibrium NWChem energy for reference
        nw_equil_file = os.path.join(args.nwchem_dir, f'water_equil.out')
        E_equil = parse_nwchem_energy(nw_equil_file)

        qi_list=[]; qj_list=[]; delta_list=[]
        HARTREE_TO_CM1 = 219474.63

        for pt in nw_points:
            nw_file = pt['file']
            E_nw    = pt['E_hartree']
            if nw_file not in job_idx:
                continue
            qi, qj = job_idx[nw_file]

            # NWChem energy relative to equilibrium (in cm-1)
            if E_equil is not None:
                V_nw = (E_nw - E_equil) * HARTREE_TO_CM1
            else:
                V_nw = E_nw * HARTREE_TO_CM1

            # MOPAC 2D GPR prediction at same point
            key = (mi, mj, round(qi,6), round(qj,6))
            if key not in mopac_data:
                continue
            V_mopac = mopac_data[key]

            delta = V_nw - V_mopac
            qi_list.append(qi); qj_list.append(qj); delta_list.append(delta)

        if len(qi_list) < 4:
            print(f"  Not enough matched points ({len(qi_list)}), skipping GPR fit")
            continue

        qi_arr = np.array(qi_list)
        qj_arr = np.array(qj_list)
        delta_arr = np.array(delta_list)

        print(f"  Matched {len(qi_arr)} points, delta range "
              f"[{delta_arr.min():.1f},{delta_arr.max():.1f}] cm-1")

        # Fit 2D delta-GPR
        gpr_delta, V_mean = fit_2d_gpr(qi_arr, qj_arr, delta_arr)
        print(f"  Fitted delta-GPR: {gpr_delta.kernel_}")

        # Save
        out_path = os.path.join(args.outdir, f'delta_gpr_2d_m{mi}m{mj}.joblib')
        joblib.dump({'gpr': gpr_delta, 'V_mean': V_mean,
                     'mode_i': mi, 'mode_j': mj}, out_path)
        print(f"  Saved {out_path}")
        results[(mi,mj)] = {'n_pts': len(qi_arr),
                             'delta_range': (delta_arr.min(), delta_arr.max())}

    print(f"\nDone. Processed {len(results)} pairs.")
    for pair, info in results.items():
        print(f"  Pair {pair}: {info['n_pts']} pts, "
              f"delta=[{info['delta_range'][0]:.1f},{info['delta_range'][1]:.1f}] cm-1")


if __name__ == '__main__':
    main()
