#!/usr/bin/env python3
"""
f34_to_morse_tanh_sop.py
Converts PBQFF f3/f4 Taylor series force constants to Morse/tanh SOP.
y_i = tanh(alpha_i * q_i), so q_i = y_i/alpha_i
Term phi * q_i^a * q_j^b * q_k^c -> phi/(alpha_i^a * alpha_j^b * alpha_k^c) * y_i^a * y_j^b * y_k^c
"""
import argparse, numpy as np
from collections import defaultdict

def load_morse(f):
    morse = {}
    for line in open(f):
        p = line.split()
        if len(p) >= 3:
            morse[int(p[0])] = (float(p[1]), float(p[2]), float(p[3]) if len(p)>3 else 0.0)
    return morse

def load_fn(f, n):
    terms = []
    for line in open(f):
        p = line.split()
        if len(p) == n+1:
            idx = tuple(int(x) for x in p[:n])
            phi = float(p[n])
            terms.append(idx + (phi,))
    return terms

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--f3",    required=True)
    p.add_argument("--f4",    required=True)
    p.add_argument("--morse", required=True)
    p.add_argument("--outfile", required=True)
    p.add_argument("--system", default="h2o")
    args = p.parse_args()

    morse = load_morse(args.morse)
    f3 = load_fn(args.f3, 3)
    f4 = load_fn(args.f4, 4)

    print(f"Loaded {len(f3)} cubic, {len(f4)} quartic terms")
    print(f"Morse alphas: {[(m,v[1]) for m,v in sorted(morse.items())]}")

    sop = []
    for terms, order in [(f3,3),(f4,4)]:
        for t in terms:
            idx = t[:order]
            phi = t[order]
            pows = defaultdict(int)
            for i in idx: pows[i] += 1
            coeff = phi
            for mode,n in pows.items():
                if mode not in morse: coeff=0; break
                alpha = morse[mode][1]
                coeff /= alpha**n
            if abs(coeff) > 0.001:
                n_unique = len(set(idx))
                sop.append((order, n_unique, dict(pows), phi, coeff, idx))

    sop.sort(key=lambda x: (x[0], x[1], -abs(x[4])))

    with open(args.outfile, "w") as f:
        f.write(f"# {args.system} Morse/tanh SOP from PBQFF f3/f4\n")
        f.write(f"# y_i = tanh(alpha_i*q_i)\n")
        for mode,(freq,alpha,shift) in sorted(morse.items()):
            f.write(f"# Mode {mode}: freq={freq:.3f} alpha={alpha:.6f}\n")
        f.write(f"#\n")
        f.write(f"# {'Order':>5} {'Modes':>5} {'Coeff (tanh basis)':>20} {'Orig phi':>15}  Powers    Indices\n")
        for order,nmodes,pows,phi,coeff,idx in sop:
            pow_str = "*".join(f"y{m}^{n}" if n>1 else f"y{m}" for m,n in sorted(pows.items()))
            idx_str = " ".join(str(i) for i in idx)
            f.write(f"  {order:5d} {nmodes:5d} {coeff:20.8f} {phi:15.8f}  {pow_str:<12}  [{idx_str}]\n")

    print(f"\nWritten {len(sop)} terms to {args.outfile}")
    print(f"\nBreakdown:")
    for order,nmodes in sorted(set((t[0],t[1]) for t in sop)):
        n = sum(1 for t in sop if t[0]==order and t[1]==nmodes)
        mx = max(abs(t[4]) for t in sop if t[0]==order and t[1]==nmodes)
        print(f"  {order}-index {nmodes}-mode: {n} terms, max|coeff|={mx:.4f} cm-1")

if __name__ == "__main__":
    main()
