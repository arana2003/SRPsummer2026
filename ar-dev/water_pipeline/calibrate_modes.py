#!/usr/bin/env python3
"""
calibrate_modes.py
Determines the correct Angstrom-per-dimensionless-q scale factor
for each normal mode by running a small MOPAC calibration.

For each mode: run MOPAC at q=0, +step, -step (step in Ang)
Then: k = 2*V/step^2, scale = sqrt(nu/k)
So that: V(q_dimless) = nu/2 * q_dimless^2 (harmonic)

Usage:
  python3 calibrate_modes.py --mol water \
      --geom pbqff_geom.npy --modes pbqff_modes.npy --freqs pbqff_freqs.npy \
      --equil-hof -57.799 --mopac /opt/mopac/mopac \
      --outfile mode_scales.npy
"""

import argparse, os, re, subprocess, numpy as np

KCAL_TO_CM1 = 349.7551

def run_mopac(geom, labels, tag, equil_hof, mopac, shm):
    mop = f"1SCF AUX(PRECISION=14 COMP) PM7 CHARGE=0\n{tag}\n{tag}\n"
    mop += "\n".join(
        f"  {lbl}  {xyz[0]:.8f} 1  {xyz[1]:.8f} 1  {xyz[2]:.8f} 1"
        for lbl, xyz in zip(labels, geom)) + "\n"
    mop_path = os.path.join(shm, tag + ".mop")
    aux_path = mop_path.replace(".mop", ".aux")
    open(mop_path, "w").write(mop)
    subprocess.run([mopac, mop_path], capture_output=True, timeout=120)
    if not os.path.exists(aux_path):
        return None
    txt = open(aux_path).read()
    m = re.search(r"HEAT_OF_FORMATION:KCAL/MOL=([\d.\-+DE]+)", txt)
    for ext in [".mop",".aux",".out",".arc",".den",".res"]:
        f = mop_path.replace(".mop", ext)
        if os.path.exists(f): os.remove(f)
    if not m: return None
    hof = float(m.group(1).replace("D","E"))
    return (hof - equil_hof) * KCAL_TO_CM1


def calibrate_mode(geom, labels, vec, freq, equil_hof, mopac, shm,
                   step=0.02, mode_idx=0):
    """
    Run MOPAC at +step and -step along mode vec (in Ang).
    Returns scale factor: step_ang = q_dimless * scale
    """
    V_pos = run_mopac(geom + step*vec, labels,
                      f"calib_m{mode_idx}_pos", equil_hof, mopac, shm)
    V_neg = run_mopac(geom - step*vec, labels,
                      f"calib_m{mode_idx}_neg", equil_hof, mopac, shm)
    if V_pos is None or V_neg is None:
        raise RuntimeError(f"Mode {mode_idx} calibration MOPAC failed")

    # Average V from both sides (removes linear term)
    V_avg = (V_pos + V_neg) / 2.0
    V_avg = max(V_avg, 0.1)  # guard against numerical noise

    # k = 2*V/step^2, scale = sqrt(nu/k)
    k = 2.0 * V_avg / step**2
    scale = np.sqrt(freq / k)

    print(f"  Mode {mode_idx+1} (freq={freq:.1f}): "
          f"V_avg={V_avg:.4f} cm-1 at step={step} Ang, "
          f"k={k:.0f} cm-1/Ang^2, scale={scale:.5f} Ang/q")
    return scale


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--mol",        default="water")
    p.add_argument("--geom",       required=True)
    p.add_argument("--modes",      required=True)
    p.add_argument("--freqs",      required=True)
    p.add_argument("--labels",     required=True)
    p.add_argument("--equil-hof",  type=float, required=True)
    p.add_argument("--mopac",      default="/opt/mopac/mopac")
    p.add_argument("--shm",        default="/dev/shm")
    p.add_argument("--outfile",    required=True)
    p.add_argument("--step",       type=float, default=0.02,
                   help="Calibration step size in Angstrom")
    args = p.parse_args()

    geom   = np.load(args.geom)
    modes  = np.load(args.modes)   # shape (n_modes, 3*n_atoms)
    freqs  = np.load(args.freqs)
    labels = list(np.load(args.labels))
    n_modes = len(freqs)

    print(f"Calibrating {n_modes} modes for {args.mol}...")
    scales = []
    for i in range(n_modes):
        vec = modes[i].reshape(-1, 3)
        scale = calibrate_mode(geom, labels, vec, freqs[i],
                               args.equil_hof, args.mopac, args.shm,
                               step=args.step, mode_idx=i)
        scales.append(scale)

    scales = np.array(scales)
    np.save(args.outfile, scales)
    print(f"\nSaved scales to {args.outfile}")
    print(f"Scales (Ang/q): {scales}")


if __name__ == "__main__":
    main()
