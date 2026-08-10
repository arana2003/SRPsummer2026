#!/usr/bin/env python3
"""
read_kaiwan_outputs.py
Reads output files from Kaiwan's pipeline:
  - nwc_nmodes_{nm}.dat  : NWChem frequencies + mode vectors (from nwchem_fc.py)
  - pbqff/pbqff.out      : PBQFF geometry and mode vectors
  - nwcf2{nm}.dat        : Matched harmonic force constants (from translation.py)

Usage:
  python3 read_kaiwan_outputs.py --nm h2o --workdir /path/to/kaiwan/output
"""

import argparse, os, numpy as np

ATOMIC_MASSES = {
    'H': 1.007825, 'C': 12.000000, 'N': 14.003074,
    'O': 15.994910, 'F': 18.998403, 'S': 31.972071,
}

def read_nwc_nmodes(nm, workdir='.'):
    """
    Read NWChem frequencies and mode vectors from nwc_nmodes_{nm}.dat
    Format: freq  vec[0]  vec[1]  ...  vec[3N-1]  (one line per mode)
    """
    path = os.path.join(workdir, 'nwchem', f'nwc_nmodes_{nm}.dat')
    if not os.path.exists(path):
        raise FileNotFoundError(f"Not found: {path}")
    freqs = []; vecs = []
    for line in open(path):
        vals = list(map(float, line.split()))
        if len(vals) < 2: continue
        freqs.append(vals[0])
        vecs.append(vals[1:])
    freqs = np.array(freqs)
    vecs  = np.array(vecs)  # shape (n_modes, 3*n_atoms)
    print(f"NWChem modes: {len(freqs)} modes, vec shape {vecs.shape}")
    return freqs, vecs


def read_pbqff_geometry(workdir='.'):
    """
    Read PBQFF optimized geometry from pbqff/pbqff.out
    PBQFF outputs geometry in BOHR -- converts to Angstrom.
    Returns (labels, coords_ang)
    """
    BOHR_TO_ANG = 0.529177
    path = os.path.join(workdir, 'pbqff', 'pbqff.out')
    if not os.path.exists(path):
        raise FileNotFoundError(f"Not found: {path}")

    labels = []; coords = []
    with open(path) as f:
        for line in f:
            if 'Geometry:' in line:
                break
        for line in f:
            if not line.strip():
                break
            vals = line.split()
            if len(vals) >= 4:
                labels.append(vals[0])
                # PBQFF geometry is in Bohr, convert to Angstrom
                coords.append([float(vals[1])*BOHR_TO_ANG,
                               float(vals[2])*BOHR_TO_ANG,
                               float(vals[3])*BOHR_TO_ANG])

    coords = np.array(coords)
    print(f"PBQFF geometry: {len(labels)} atoms {labels} (converted Bohr->Ang)")
    return labels, coords


def read_pbqff_nmodes(nm, workdir='.'):
    path = os.path.join(workdir, 'pbqff', f'pbqff2_nmodes_{nm}.dat')
    print(f"Reading PBQFF modes: {path}")
    return read_nmodes_file(path)


def read_nwchem_geometry(nm, workdir='.'):
    """
    Read NWChem optimized geometry from nwchem/{nm}.out
    Returns coords_ang
    """
    path = os.path.join(workdir, 'nwchem', f'nwc_{nm}.out')
    if not os.path.exists(path):
        # Try alternate locations
        for alt in [f'{nm}_opt.out', f'{nm}.out']:
            alt_path = os.path.join(workdir, alt)
            if os.path.exists(alt_path):
                path = alt_path
                break

    coords = []
    with open(path) as f:
        for line in f:
            if 'Output coordinates in angstroms' in line:
                # Skip header lines
                for _ in range(3): next(f)
                break
        for line in f:
            if not line.strip():
                break
            vals = line.split()
            if len(vals) >= 6 and vals[0].isdigit():
                coords.append([float(vals[3]), float(vals[4]), float(vals[5])])

    coords = np.array(coords)
    print(f"NWChem geometry: {len(coords)} atoms")
    return coords


def save_as_npy(nm, workdir, outdir):
    """Extract all data and save as .npy files for worker_v3.py"""
    os.makedirs(outdir, exist_ok=True)

    # MOPAC path: PBQFF geometry + PBQFF modes
    labels, pbqff_geom = read_pbqff_geometry(workdir)
    masses = np.array([ATOMIC_MASSES.get(l, 12.0) for l in labels])

    # NWChem path: NWChem geometry + NWChem modes
    nwc_freqs, nwc_vecs = read_nwc_nmodes(nm, workdir)
    nwc_geom = read_nwchem_geometry(nm, workdir)

    # Save MOPAC/PBQFF path files
    np.save(f'{outdir}/{nm}_pbqff_geometry.npy', pbqff_geom)
    np.save(f'{outdir}/{nm}_labels.npy', np.array(labels))
    np.save(f'{outdir}/{nm}_masses.npy', masses)

    # Save NWChem path files
    np.save(f'{outdir}/{nm}_nwchem_geometry.npy', nwc_geom)
    np.save(f'{outdir}/{nm}_nwchem_frequencies.npy', nwc_freqs)
    np.save(f'{outdir}/{nm}_nwchem_modes.npy', nwc_vecs)

    print(f"\nSaved all .npy files to {outdir}/")
    print(f"  PBQFF geometry: {pbqff_geom.shape}")
    print(f"  NWChem modes:   {nwc_vecs.shape}")
    print(f"  NWChem freqs:   {nwc_freqs}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--nm',      required=True, help='System name e.g. h2o')
    p.add_argument('--workdir', default='.', help='Kaiwan pipeline working directory')
    p.add_argument('--outdir',  default='molecule_data')
    args = p.parse_args()
    save_as_npy(args.nm, args.workdir, args.outdir)


if __name__ == '__main__':
    main()
