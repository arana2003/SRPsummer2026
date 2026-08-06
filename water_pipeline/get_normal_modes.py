#!/usr/bin/env python3
"""
get_normal_modes.py
Runs a MOPAC FORCE calculation at the equilibrium geometry to get
normal mode vectors in MOPAC's convention. These are the correct
vectors to use for displacement calculations in worker_v3.py.

Usage:
  python3 get_normal_modes.py --geom geometry.xyz --mol water \
      --mopac /opt/mopac/mopac --outdir /path/to/output

Output:
  normal_modes.npy   -- shape (n_modes, 3*n_atoms), mass-weighted
  frequencies.npy    -- shape (n_modes,) in cm-1
  geometry.npy       -- shape (n_atoms, 3) in Angstrom
  atom_labels.npy    -- shape (n_atoms,) string array
"""

import argparse, os, re, subprocess, tempfile
import numpy as np

MOPAC_FORCE_TEMPLATE = """\
FORCE AUX(PRECISION=14) PM7 CHARGE=0
{title}
normal mode calculation
{geometry}
"""

ATOMIC_MASSES = {
    'H': 1.007825, 'C': 12.000000, 'N': 14.003074,
    'O': 15.994910, 'F': 18.998403, 'S': 31.972071,
    'Cl': 34.968853, 'Br': 78.918338,
}

def read_xyz(xyz_file):
    """Read .xyz file -> (labels, coords_ang)"""
    lines = open(xyz_file).readlines()
    n = int(lines[0].strip())
    labels = []; coords = []
    for line in lines[2:2+n]:
        parts = line.split()
        labels.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return labels, np.array(coords)


def write_mopac_input(path, title, labels, coords):
    geom_lines = [
        f"  {lbl:4s}  {xyz[0]:14.8f} 1  {xyz[1]:14.8f} 1  {xyz[2]:14.8f} 1"
        for lbl, xyz in zip(labels, coords)
    ]
    with open(path, 'w') as f:
        f.write(MOPAC_FORCE_TEMPLATE.format(
            title=title, geometry="\n".join(geom_lines)))


def parse_mopac_force_output(aux_path, n_atoms, freq_cutoff=200.0):
    """Parse frequencies and normal mode vectors from MOPAC FORCE aux file."""
    if not os.path.exists(aux_path):
        raise FileNotFoundError(f"aux file not found: {aux_path}")

    txt = open(aux_path).read()

    # Parse frequencies
    freq_m = re.search(r'VIB\._FREQ:CM\(-1\)\[\d+\]=([\s\d.\-+DE\n]+?)(?:\n [A-Z]|\Z)', txt)
    if not freq_m:
        raise ValueError("Frequencies not found in aux file")
    freqs_all = [float(x.replace('D','E')) for x in freq_m.group(1).split()]

    # Parse normal mode vectors (FORCE_CONSTANTS or NORMAL_MODES)
    vec_m = re.search(r'NORMAL_MODES\[\d+\]=([\s\d.\-+DE\n]+?)(?:\n #|\n [A-Z]|\Z)', txt)
    if not vec_m:
        raise ValueError("Normal modes not found in aux file")
    vecs_flat = [float(x.replace('D','E')) for x in vec_m.group(1).split()]

    n_coords = 3 * n_atoms
    n_modes_all = len(freqs_all)
    vecs = np.array(vecs_flat).reshape(n_modes_all, n_coords)

    # Remove translation/rotation (frequencies near zero)
    real_mask = np.array(freqs_all) > freq_cutoff
    freqs = np.array(freqs_all)[real_mask]
    vecs  = vecs[real_mask]

    print(f"Found {len(freqs)} real frequencies (removed {real_mask.size - real_mask.sum()} near-zero)")
    return freqs, vecs


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--geom',   required=True, help='.xyz file with equilibrium geometry')
    p.add_argument('--mol',    default='molecule')
    p.add_argument('--mopac',  default='/opt/mopac/mopac')
    p.add_argument('--outdir', default='.')
    p.add_argument('--charge', type=int, default=0)
    p.add_argument('--freq-cutoff', type=float, default=200.0,
                   help='Min frequency cm-1 to keep (removes translations/rotations)')
    p.add_argument('--shm',    default='/dev/shm')
    args = p.parse_args()

    # Read geometry
    labels, coords = read_xyz(args.geom)
    n_atoms = len(labels)
    masses  = np.array([ATOMIC_MASSES.get(l, 12.0) for l in labels])
    print(f"Molecule: {args.mol}, {n_atoms} atoms: {labels}")

    # Step 1: Run MOPAC geometry optimization
    opt_path = os.path.join(args.shm, f"{args.mol}_opt.mop")
    opt_aux  = opt_path.replace('.mop', '.aux')

    # Write optimization input (all coords free to optimize)
    geom_lines_opt = [
        f"  {lbl:4s}  {xyz[0]:14.8f} 1  {xyz[1]:14.8f} 1  {xyz[2]:14.8f} 1"
        for lbl, xyz in zip(labels, coords)
    ]
    with open(opt_path, 'w') as f:
        f.write(f"PM7 CHARGE={args.charge} AUX(PRECISION=14)\n{args.mol}_opt\n{args.mol}_opt\n")
        f.write("\n".join(geom_lines_opt) + "\n")

    print("Running MOPAC geometry optimization...")
    subprocess.run([args.mopac, opt_path], capture_output=True, timeout=600)

    # Extract optimized geometry - take LAST ATOM_X_UPDATED block
    import re
    opt_txt = open(opt_aux).read()
    all_geoms = re.findall(
        r"ATOM_X_UPDATED:ANGSTROMS\[\d+\]=([\s\d.\-+DE\n]+?)(?:\n [A-Z]|\Z)",
        opt_txt)
    if all_geoms:
        opt_vals = [float(x.replace("D","E")) for x in all_geoms[-1].split()]
        coords = np.array(opt_vals).reshape(-1, 3)
        print(f"Optimized geometry ({len(all_geoms)} steps):")
        for lbl,xyz in zip(labels,coords):
            print(f"  {lbl}  {xyz[0]:.6f}  {xyz[1]:.6f}  {xyz[2]:.6f}")
    else:
        print("WARNING: Could not find optimized geometry")

    # Get final HOF
    all_hofs = re.findall(r"HEAT_OF_FORM_UPDATED:KCAL/MOL=([\d.\-+DE]+)", opt_txt)
    if all_hofs:
        hof_equil = float(all_hofs[-1].replace("D","E"))
        print(f"Equilibrium HOF: {hof_equil:.6f} kcal/mol")
        os.makedirs(args.outdir, exist_ok=True)
        np.save(os.path.join(args.outdir, f"{args.mol}_equil_hof.npy"), [hof_equil])

    # Step 2: Run MOPAC FORCE at optimized geometry
    mop_path = os.path.join(args.shm, f"{args.mol}_force.mop")
    aux_path = mop_path.replace('.mop', '.aux')
    write_mopac_input(mop_path, args.mol, labels, coords)

    print("Running MOPAC FORCE calculation...")
    r = subprocess.run([args.mopac, mop_path], capture_output=True, timeout=300)
    if not os.path.exists(aux_path):
        raise RuntimeError(f"MOPAC failed: {r.stderr.decode()[:200]}")
    print("MOPAC FORCE done.")

    # Parse output
    freqs, vecs = parse_mopac_force_output(aux_path, n_atoms, args.freq_cutoff)
    print(f"Frequencies (cm-1): {freqs}")
    print(f"Mode vectors shape: {vecs.shape}")

    # Verify dimensionless convention:
    # For mode i, displacement of q=1 should give energy increase of omega/2
    # omega = freq_i / HARTREE_TO_CM1
    # In MOPAC convention, vec is mass-weighted so:
    # step_bohr = q / sqrt(omega) / mw_norm
    # where mw_norm = sqrt(sum(m_i * L_i^2))
    HARTREE_TO_CM1 = 219474.63
    AMU_TO_AU = 1822.888486
    print("\nVerification (q=1 should give ~omega/2 energy increase):")
    for i, (freq, vec) in enumerate(zip(freqs[:3], vecs[:3])):
        omega = freq / HARTREE_TO_CM1
        masses_au = masses * AMU_TO_AU
        mw_norm = np.sqrt(np.sum(masses_au[:, None] * vec.reshape(n_atoms, 3)**2))
        step_bohr = 1.0 / np.sqrt(omega) / mw_norm
        print(f"  Mode {i+1} (freq={freq:.1f}): mw_norm={mw_norm:.4f}, "
              f"q=1 -> {step_bohr:.4f} Bohr displacement")

    # Save
    os.makedirs(args.outdir, exist_ok=True)
    np.save(os.path.join(args.outdir, f'{args.mol}_geometry.npy'), coords)
    np.save(os.path.join(args.outdir, f'{args.mol}_labels.npy'), np.array(labels))
    np.save(os.path.join(args.outdir, f'{args.mol}_frequencies.npy'), freqs)
    np.save(os.path.join(args.outdir, f'{args.mol}_normal_modes.npy'), vecs)
    print(f"\nSaved to {args.outdir}/")

    # Cleanup
    for ext in ['.mop','.aux','.out','.arc']:
        f = mop_path.replace('.mop', ext)
        if os.path.exists(f): os.remove(f)


if __name__ == '__main__':
    main()
