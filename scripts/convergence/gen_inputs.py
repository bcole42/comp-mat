"""Generate every input file for the Si convergence tests.

Creates:
  conv_ecut/si8_ec{E}.in    8-atom rattled cell, scan of ecutwfc   (fast)
  conv_k/si64_k{N}.in       64-atom rattled cell, scan of k-mesh   (slower)
"""
import os
from search_common import *

os.makedirs("conv_ecut", exist_ok=True)
os.makedirs("conv_k", exist_ok=True)

# --- cutoff scan: 8-atom cell ---
cell, pos = build_si(n=1)
for ec in ECUT_LIST:
    write_pw_input(f"conv_ecut/si8_ec{ec}.in", f"si8_ec{ec}", cell, pos, ec, K_FOR_ECUT_SCAN)

# --- k-point scan: 64-atom cell ---
cell, pos = build_si(n=2)
for k in K_LIST:
    write_pw_input(f"conv_k/si64_k{k}.in", f"si64_k{k}", cell, pos, ECUT_FOR_K_SCAN, k)

print("Wrote", len(ECUT_LIST), "cutoff inputs and", len(K_LIST), "k-point inputs.")
