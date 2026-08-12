#!/usr/bin/env python3
"""
generate_job_list.py
Generates a flat CSV job list for all 1D and 2D MOPAC calculations.
Each row is one MOPAC job: (job_id, mode_i, mode_j, qi, qj)
Workers read this file and process their strided slice.

Usage:
  python3 generate_job_list.py --mol water --n-modes 3 --type 1d --outfile jobs_1d.csv
  python3 generate_job_list.py --mol methanol --n-modes 12 --type 2d --outfile jobs_2d.csv
  python3 generate_job_list.py --mol methanol --n-modes 12 --type all --outfile jobs_all.csv
"""

import argparse
import csv
import numpy as np
from itertools import combinations

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--mol',      default='water')
    p.add_argument('--n-modes',  type=int, required=True)
    p.add_argument('--type',     choices=['1d','2d','all'], default='all')
    p.add_argument('--n-1d',     type=int, default=55,  help='Points per 1D cut')
    p.add_argument('--n-2d',     type=int, default=7,   help='Grid size per 2D cut (n x n)')
    p.add_argument('--q-max-1d', type=float, default=1.7, help='Max |q| for 1D')
    p.add_argument('--q-max-2d', type=float, default=1.7, help='Max |q| for 2D')
    p.add_argument('--outfile',  required=True)
    args = p.parse_args()

    jobs = []
    job_id = 0

    # 1D cuts: displace along one mode at a time (qj=0)
    if args.type in ['1d', 'all']:
        q_1d = np.linspace(-args.q_max_1d, args.q_max_1d, args.n_1d)
        for mode in range(1, args.n_modes + 1):
            for qi in q_1d:
                jobs.append({
                    'job_id':  job_id,
                    'mol':     args.mol,
                    'mode_i':  mode,
                    'mode_j':  mode,
                    'qi':      round(float(qi), 8),
                    'qj':      0.0,
                    'type':    '1d',
                })
                job_id += 1
        print(f"1D jobs: {args.n_modes} modes x {args.n_1d} pts = {args.n_modes * args.n_1d}")

    # 2D cuts: displace along two modes simultaneously
    if args.type in ['2d', 'all']:
        q_2d = np.linspace(-args.q_max_2d, args.q_max_2d, args.n_2d)
        pairs = list(combinations(range(1, args.n_modes + 1), 2))
        for mi, mj in pairs:
            for qi in q_2d:
                for qj in q_2d:
                    jobs.append({
                        'job_id':  job_id,
                        'mol':     args.mol,
                        'mode_i':  mi,
                        'mode_j':  mj,
                        'qi':      round(float(qi), 8),
                        'qj':      round(float(qj), 8),
                        'type':    '2d',
                    })
                    job_id += 1
        n_pairs = len(pairs)
        n_2d_pts = args.n_2d ** 2
        print(f"2D jobs: {n_pairs} pairs x {n_2d_pts} pts = {n_pairs * n_2d_pts}")

    # Write CSV
    fields = ['job_id','mol','mode_i','mode_j','qi','qj','type']
    with open(args.outfile, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(jobs)

    print(f"Total: {len(jobs)} jobs -> {args.outfile}")

if __name__ == '__main__':
    main()
