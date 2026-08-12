#!/usr/bin/env python3
"""
parse_sample_inp.py
Parses Kaiwan's sample.inp file to extract molecule info
needed by the GPR pipeline.

Usage:
  python3 parse_sample_inp.py --inp sample.inp --vars system.vars
"""

import argparse, re, os, numpy as np

def parse_system_vars(vars_file):
    """Parse system.vars shell export file."""
    env = {}
    for line in open(vars_file):
        m = re.match(r'export\s+(\w+)=["\']?([^"\']+)["\']?', line.strip())
        if m:
            env[m.group(1)] = m.group(2).strip()
    return env

def parse_sample_inp(inp_file):
    """Parse sample.inp into sections separated by --"""
    content = open(inp_file).read()
    sections = [s.strip() for s in content.split('--')]

    result = {}

    # Section 0: system_name and geometry
    sec0 = sections[0]
    name_m = re.search(r'system_name\s*=\s*(\S+)', sec0)
    result['system_name'] = name_m.group(1) if name_m else 'unknown'

    # Parse geometry block (atom lines)
    geom_lines = []
    labels = []
    for line in sec0.splitlines():
        parts = line.split()
        if len(parts) == 4 and parts[0].isalpha():
            labels.append(parts[0])
            geom_lines.append([float(parts[1]), float(parts[2]), float(parts[3])])
    result['labels']   = labels
    result['geometry'] = np.array(geom_lines) if geom_lines else None

    # Section 1: PBQFF settings
    if len(sections) > 1:
        pbqff = sections[1]
        result['charge']    = int(re.search(r'charge\s*=\s*(\d+)', pbqff).group(1)) \
                              if re.search(r'charge\s*=\s*(\d+)', pbqff) else 0
        result['program']   = re.search(r'program\s*=\s*["\']?(\w+)', pbqff).group(1) \
                              if re.search(r'program\s*=\s*["\']?(\w+)', pbqff) else 'mopac'
        result['optimize']  = 'true' in pbqff.lower()
        result['coord_type']= re.search(r'coord_type\s*=\s*["\']?(\w+)', pbqff).group(1) \
                              if re.search(r'coord_type\s*=\s*["\']?(\w+)', pbqff) else 'normal'
        result['step_size'] = float(re.search(r'step_size\s*=\s*([\d.]+)', pbqff).group(1)) \
                              if re.search(r'step_size\s*=\s*([\d.]+)', pbqff) else 0.005

    # Section 4: MLCP settings (last section)
    if len(sections) > 4:
        result['mlcp_input'] = sections[4]

    return result

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--inp',  required=True, help='sample.inp file')
    p.add_argument('--vars', default=None,  help='system.vars file')
    p.add_argument('--outdir', default='.', help='Output directory for .npy files')
    args = p.parse_args()

    # Parse system.vars
    env = {}
    if args.vars and os.path.exists(args.vars):
        env = parse_system_vars(args.vars)
        print("System vars:")
        for k,v in env.items():
            print(f"  {k} = {v}")

    # Parse sample.inp
    info = parse_sample_inp(args.inp)
    print(f"\nSystem name: {info['system_name']}")
    print(f"Labels: {info['labels']}")
    print(f"Geometry (Angstrom):\n{info['geometry']}")
    print(f"Program: {info['program']}")
    print(f"Charge: {info['charge']}")
    print(f"Coord type: {info['coord_type']}")
    print(f"Step size: {info['step_size']}")

    # Save geometry and labels
    if info['geometry'] is not None:
        os.makedirs(args.outdir, exist_ok=True)
        name = info['system_name']
        np.save(f"{args.outdir}/{name}_geometry.npy", info['geometry'])
        np.save(f"{args.outdir}/{name}_labels.npy", np.array(info['labels']))
        print(f"\nSaved {name}_geometry.npy and {name}_labels.npy to {args.outdir}/")

if __name__ == '__main__':
    main()
