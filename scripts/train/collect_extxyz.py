"""Stage 3: turn the QE outputs into train.extxyz / valid.extxyz for MACE fine-tuning.

  python3 03_collect_extxyz.py              # default: align energy offset to the foundation model
  python3 03_collect_extxyz.py --no-align   # keep raw QE total energies

Keys written (pass them explicitly to mace_run_train):
  --energy_key="REF_energy" --forces_key="REF_forces" --stress_key="REF_stress"

Conventions:
  energy  eV, forces eV/A, stress eV/A^3 in the ASE sign convention
  (QE prints pressure-like stress, ASE uses the opposite sign, so stress_ASE = -stress_QE).

Energy alignment: QE and the VASP-PAW data behind MACE-MP-0 use different pseudopotentials,
so absolute total energies differ by a large constant per atom. With --align (default) a single
constant shift per atom is subtracted so the mean DFT energy matches the foundation model's
energies on the same structures (fit on training frames only). Forces and stress are untouched,
and these are what matter most for phonons. This is a pragmatic fix for one element; check how
your fine-tuning recipe handles isolated-atom energies (E0s).
"""
import json
import sys
import numpy as np
from train_common import *

align = "--no-align" not in sys.argv
manifest = json.load(open("manifest.json"))

results, bad = [], []
for m in manifest:
    try:
        results.append(parse_pw_output(f"qe/{m['name']}.out"))
    except Exception as err:
        bad.append(f"{m['name']}: {err}")
if bad:
    sys.exit("Problems with:\n  " + "\n  ".join(bad))

# per-atom offset between DFT and foundation model, from training frames only
d = np.array([(r["energy"] - m["mace_energy"]) / m["n_atoms"]
              for r, m in zip(results, manifest) if m["split"] == "train"])
delta = float(d.mean()) if align else 0.0
print(f"DFT - foundation-model energy offset: {d.mean():.4f} eV/atom "
      f"(spread {d.std() * 1000:.1f} meV/atom over training frames)")
print("Applying shift" if align else "Not applying shift (--no-align)")

frames = {"train": [], "valid": []}
for r, m in zip(results, manifest):
    n = m["n_atoms"]
    cell = np.array(m["cell"])
    vol = abs(np.linalg.det(cell))
    frames[m["split"]].append(dict(
        cell=cell, pos=np.array(m["pos"]), forces=r["forces"],
        energy=r["energy"] - n * delta,
        stress=-r["stress_qe"] * KBAR_TO_EV_A3,
        info=dict(config_type=m["split"], vol_pct=m["vol_pct"], name=m["name"],
                  dft_energy_raw=f"{r['energy']:.8f}", volume_per_atom=f"{vol / n:.5f}")))

for split, fr in frames.items():
    write_extxyz(f"{split}.extxyz", fr)
    fmax = max(np.abs(f["forces"]).max() for f in fr)
    frms = np.sqrt(np.mean(np.concatenate([f["forces"].ravel() for f in fr]) ** 2))
    print(f"{split}.extxyz: {len(fr)} frames, force RMS {frms * 1000:.0f} meV/A, max {fmax:.3f} eV/A")

print("\nSign check: the -5% volume frames should have NEGATIVE stress diagonals (compressed),")
print("the +5% volume frames POSITIVE ones (stretched), in the ASE convention.")
seen = set()
for fr in frames["train"]:
    v = fr["info"]["vol_pct"]
    if v in (-5.0, 5.0) and v not in seen:
        seen.add(v)
        print(f"  vol {v:+.1f}%  mean stress diag = "
              f"{np.trace(fr['stress']) / 3 / KBAR_TO_EV_A3:+.1f} kbar")
