#!/usr/bin/env python3
import argparse, csv, fcntl, os, re, subprocess, time
import numpy as np

AMU_TO_AU      = 1822.888486
BOHR_TO_ANG    = 0.529177
HARTREE_TO_CM1 = 219474.63
KCAL_TO_CM1    = 349.7551
KCAL_TO_AU     = 1.0 / 627.5094740631
ANG_TO_BOHR    = 1.0 / BOHR_TO_ANG

CSV_FIELDS = ['job_id','mol','mode_i','mode_j','qi','qj',
              'hof_kcal','V_cm1','dV_dqi','dV_dqj','status']

MOPAC_TEMPLATE = ("1SCF GRADIENTS AUX(PRECISION=14 COMP) PM7 CHARGE=0\n"
                  "{title}\nworker_v3\n{geometry}\n")

def load_molecule(mol_name):
    if mol_name == 'water':
        d=os.environ.get('MOLDATA_DIR','/global/homes/a/alexisr/my_container_build/water_widerange_test/wide_range_mopac_inputs')+'/' 
        geom    = np.load(d+'geometry.npy')
        labels  = list(np.load(d+'atom_labels.npy'))
        modes   = np.load(d+'normal_modes.npy')
        freqs   = np.load(d+'frequencies.npy')
        masses  = np.array([15.9994, 1.00794, 1.00794])
        return {
            'labels':    labels,
            'masses':    masses,
            'geom_ang':  geom,
            'mode_vecs': {i+1: modes[i] for i in range(len(freqs))},
            'freqs':     {i+1: freqs[i] for i in range(len(freqs))},
        }
    else:
        raise NotImplementedError(f"Add {mol_name} to load_molecule()")

def precompute(mol, mi, mj):
    masses_au = mol['masses'] * AMU_TO_AU
    mass_vec  = np.repeat(np.sqrt(masses_au), 3)
    def mw_norm(vec):
        return np.sqrt(np.sum(masses_au[:, None] * vec.reshape(-1, 3)**2))
    return {
        'mass_vec':  mass_vec,
        'freq_i':    mol['freqs'][mi],
        'freq_j':    mol['freqs'][mj],
        'omega_i':   mol['freqs'][mi] / HARTREE_TO_CM1,
        'omega_j':   mol['freqs'][mj] / HARTREE_TO_CM1,
        'vec_i':     mol['mode_vecs'][mi],
        'vec_j':     mol['mode_vecs'][mj],
        'n_atoms':   len(mol['labels']),
    }

def q_to_bohr(q, freq_cm1, eigenvector_per_atom, atom_masses_amu):
    omega_hartree = freq_cm1 / HARTREE_TO_CM1
    masses_au = np.asarray(atom_masses_amu) * AMU_TO_AU
    mw_norm = np.sqrt(np.sum(masses_au[:, None] * eigenvector_per_atom**2))
    return (q / np.sqrt(omega_hartree)) / mw_norm


def make_geometry(mol, pre, qi, qj):
    n = pre['n_atoms']
    vec_i = pre['vec_i'].reshape(n, 3)
    vec_j = pre['vec_j'].reshape(n, 3)
    step_i = q_to_bohr(qi, pre['freq_i'], vec_i, mol['masses'])
    step_j = q_to_bohr(qj, pre['freq_j'], vec_j, mol['masses'])
    geom = (mol['geom_ang']
            + step_i * BOHR_TO_ANG * vec_i
            + step_j * BOHR_TO_ANG * vec_j)
    lines = [f"  {lbl:4s}  {xyz[0]:14.8f} 1  {xyz[1]:14.8f} 1  {xyz[2]:14.8f} 1"
             for lbl, xyz in zip(mol['labels'], geom)]
    return "\n".join(lines)

def parse_aux(path, pre, equil_hof):
    if not os.path.exists(path):
        return None, None, None, None, "aux file missing"
    txt = open(path).read()
    if 'UNABLE TO ACHIEVE SELF-CONSISTENT' in txt or 'SCF FAILED' in txt:
        return None, None, None, None, "SCF convergence failure"
    if 'HEAT_OF_FORMATION' not in txt:
        return None, None, None, None, "abnormal termination"
    m = re.search(r'HEAT_OF_FORMATION:KCAL/MOL=([\d.\-+DE]+)', txt)
    if not m:
        return None, None, None, None, "HOF not found"
    hof   = float(m.group(1).replace('D', 'E'))
    V_cm1 = (hof - equil_hof) * KCAL_TO_CM1
    dV_dqi = dV_dqj = None
    gm = re.search(r'GRADIENTS:KCAL/MOL/ANGSTROM\[\d+\]=\n([^\n]+)', txt)
    if gm:
        try:
            vals   = [float(x.replace('D','E')) for x in gm.group(1).split()]
            g_au   = np.array(vals) * KCAL_TO_AU * ANG_TO_BOHR
            g_mw   = g_au / pre['mass_vec']
            dV_dqi = np.dot(g_mw, pre['vec_i']) / np.sqrt(pre['omega_i']) * HARTREE_TO_CM1
            dV_dqj = np.dot(g_mw, pre['vec_j']) / np.sqrt(pre['omega_j']) * HARTREE_TO_CM1
        except Exception:
            pass
    return hof, V_cm1, dV_dqi, dV_dqj, None

def cleanup_shm(shm_dir, tag):
    for ext in ['.mop','.aux','.out','.arc','.den','.res','.pdb','.end']:
        f = os.path.join(shm_dir, tag + ext)
        if os.path.exists(f):
            try: os.remove(f)
            except OSError: pass

def append_row(csv_path, row):
    for attempt in range(30):
        try:
            with open(csv_path, 'a', newline='') as f:
                fcntl.flock(f, fcntl.LOCK_EX)
                empty = (os.path.getsize(csv_path) == 0)
                w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
                if empty: w.writeheader()
                w.writerow(row)
                fcntl.flock(f, fcntl.LOCK_UN)
            return
        except (BlockingIOError, OSError):
            time.sleep(0.01 * (attempt + 1))
    raise RuntimeError(f"Cannot write to {csv_path}")

def append_retry(retry_path, qi, qj, reason):
    with open(retry_path, 'a') as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        f.write(f"{qi:.8f} {qj:.8f}  # {reason}\n")
        fcntl.flock(f, fcntl.LOCK_UN)

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--worker-id',     type=int, required=True)
    p.add_argument('--n-workers',     type=int, default=128)
    p.add_argument('--node-id',       type=int, default=0)
    p.add_argument('--jobs-per-node', type=int, default=32000)
    p.add_argument('--mol',           default='water')
    p.add_argument('--mode-i',        type=int, required=True)
    p.add_argument('--mode-j',        type=int, required=True)
    p.add_argument('--pointsfile',    required=True)
    p.add_argument('--outcsv',        required=True)
    p.add_argument('--equil-hof',     type=float, default=0.0)
    p.add_argument('--mopac',         default='/opt/mopac/mopac')
    p.add_argument('--shm',           default='/dev/shm')
    p.add_argument('--time-limit',    type=float, default=None)
    args = p.parse_args()

    t_start = time.time()
    mi, mj  = args.mode_i, args.mode_j
    retry_path = args.outcsv.replace('.csv',
                     f'_retry_n{args.node_id}w{args.worker_id}.txt')

    mol = load_molecule(args.mol)
    pre = precompute(mol, mi, mj)

    with open(args.pointsfile) as f:
        all_points = [tuple(map(float, ln.split()[:2]))
                      for ln in f if ln.strip() and not ln.startswith('#')]

    node_start  = args.node_id * args.jobs_per_node
    node_end    = min(node_start + args.jobs_per_node, len(all_points))
    node_points = all_points[node_start:node_end]
    my_points   = node_points[args.worker_id::args.n_workers]

    print(f"Worker n{args.node_id}w{args.worker_id}: "
          f"{len(my_points)} pts assigned")

    n_done = n_fail = n_skipped = 0

    for local_idx, (qi, qj) in enumerate(my_points):
        if args.time_limit and (time.time() - t_start) > args.time_limit:
            for qi2, qj2 in my_points[local_idx:]:
                append_retry(retry_path, qi2, qj2, "time_limit")
            n_skipped = len(my_points) - local_idx
            break

        global_idx = node_start + args.worker_id + local_idx * args.n_workers
        tag     = f"n{args.node_id}w{args.worker_id}j{global_idx}_m{mi}m{mj}"
        shm_mop = os.path.join(args.shm, tag + '.mop')
        shm_aux = os.path.join(args.shm, tag + '.aux')

        try:
            geom = make_geometry(mol, pre, qi, qj)
            with open(shm_mop, 'w') as f:
                f.write(MOPAC_TEMPLATE.format(title=tag, geometry=geom))
            subprocess.run([args.mopac, shm_mop],
                           capture_output=True, timeout=300)
            hof, V_cm1, dV_dqi, dV_dqj, err = parse_aux(
                shm_aux, pre, args.equil_hof)
            if err:
                raise RuntimeError(err)
            append_row(args.outcsv, {
                'job_id': global_idx, 'mol': args.mol,
                'mode_i': mi, 'mode_j': mj, 'qi': qi, 'qj': qj,
                'hof_kcal': round(hof, 8), 'V_cm1': round(V_cm1, 4),
                'dV_dqi': round(dV_dqi, 4) if dV_dqi is not None else '',
                'dV_dqj': round(dV_dqj, 4) if dV_dqj is not None else '',
                'status': 'ok',
            })
            n_done += 1
        except subprocess.TimeoutExpired:
            append_retry(retry_path, qi, qj, "mopac_timeout")
            n_fail += 1
        except Exception as e:
            reason = str(e) or "unknown"
            append_retry(retry_path, qi, qj, reason)
            try:
                append_row(args.outcsv, {
                    'job_id': global_idx, 'mol': args.mol,
                    'mode_i': mi, 'mode_j': mj, 'qi': qi, 'qj': qj,
                    'hof_kcal':'','V_cm1':'','dV_dqi':'','dV_dqj':'',
                    'status': f'FAILED:{reason}',
                })
            except Exception:
                pass
            n_fail += 1
        finally:
            cleanup_shm(args.shm, tag)

    elapsed = time.time() - t_start
    print(f"Worker n{args.node_id}w{args.worker_id}: "
          f"{n_done} ok | {n_fail} failed | {n_skipped} skipped | "
          f"{elapsed:.1f}s")

if __name__ == '__main__':
    main()
