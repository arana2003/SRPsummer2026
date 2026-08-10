#!/usr/bin/env python3
"""
detect_symmetry.py
Automatically determines if each mode has a symmetric or asymmetric potential
by comparing V(+q) vs V(-q) at a few test displacements.

Symmetric modes  -> tanh basis (even powers)
Asymmetric modes -> Morse basis

Usage:
  python3 detect_symmetry.py --results-csv results_1d.csv --threshold 50
"""

import argparse, csv, numpy as np

def detect_symmetry(results_csv, threshold_cm1=50.0):
    """
    For each mode, compare V(+q) and V(-q) at matching displacements.
    If max |V(+q) - V(-q)| > threshold -> asymmetric (Morse)
    Otherwise -> symmetric (tanh)
    """
    rows = [r for r in csv.DictReader(open(results_csv)) if r['status']=='ok']

    # Group by mode
    modes = {}
    for r in rows:
        mi = int(r['mode_i'])
        qi = float(r['qi'])
        V  = float(r['V_cm1'])
        if mi not in modes:
            modes[mi] = {}
        modes[mi][round(qi, 6)] = V

    results = {}
    for mode, data in sorted(modes.items()):
        q_vals = sorted(data.keys())
        V_ref  = min(data.values())

        # Find matching +/- pairs
        asymmetries = []
        for q in q_vals:
            if q > 0 and round(-q, 6) in data:
                V_pos = data[q]      - V_ref
                V_neg = data[round(-q,6)] - V_ref
                asymmetries.append(abs(V_pos - V_neg))

        if not asymmetries:
            results[mode] = 'symmetric'
            continue

        max_asym = max(asymmetries)
        is_sym   = max_asym < threshold_cm1
        results[mode] = 'symmetric' if is_sym else 'asymmetric'
        print(f"Mode {mode:2d}: max asymmetry = {max_asym:8.2f} cm-1 -> {results[mode]}")

    sym_modes  = [m for m,s in results.items() if s=='symmetric']
    asym_modes = [m for m,s in results.items() if s=='asymmetric']
    print(f"\nSymmetric  (tanh basis): modes {sym_modes}")
    print(f"Asymmetric (Morse basis): modes {asym_modes}")
    return results

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--results-csv', required=True)
    p.add_argument('--threshold',   type=float, default=50.0,
                   help='Max |V(+q)-V(-q)| in cm-1 to call symmetric')
    args = p.parse_args()
    detect_symmetry(args.results_csv, args.threshold)

if __name__ == '__main__':
    main()
