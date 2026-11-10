"""Compare each run to the most converged one (highest cutoff / densest k-mesh).

Targets: max force change < ~1 meV/A, stress change ~< 1 kbar.
"""
import numpy as np
from conv_common import *


def report(label, files, values):
    res = [parse_pw_output(f) for f in files]
    ref = res[-1]
    n = len(ref["forces"])
    print(f"\n{label}   (reference = last row, {n} atoms)")
    print(f"{'value':>6} {'max|dF| meV/A':>15} {'max|dStress| kbar':>19} {'dE meV/atom':>13}")
    for v, r in zip(values, res):
        dF = np.abs(r["forces"] - ref["forces"]).max() * 1000
        dS = np.abs(r["stress"] - ref["stress"]).max()
        dE = (r["energy"] - ref["energy"]) / n * 1000
        flag = "  <-- ok" if dF < 1.0 and dS < 1.0 else ""
        print(f"{v:>6} {dF:>15.3f} {dS:>19.3f} {dE:>13.3f}{flag}")


report("ecutwfc scan (Ry), 8-atom cell",
       [f"conv_ecut/si8_ec{e}.out" for e in ECUT_LIST], ECUT_LIST)
report("k-mesh scan (NxNxN), 64-atom cell",
       [f"conv_k/si64_k{k}.out" for k in K_LIST], K_LIST)
