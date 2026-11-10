import re
import numpy as np

# Settings
A_LATTICE = 5.469            # Angstrom, PBE lattice constant of diamond Si
PSEUDO_FILE = "Si.us.pbesol.z_4.ld1.psl.v1.0.0-high.upf"
PSEUDO_DIR = "../pseudo"
MASS_SI = 28.0855
DUAL = 8                     # ecutrho = DUAL * ecutwfc (ultrasoft/PAW)
CONV_THR = "1.0d-10"         # Ry, whole cell
RATTLE_STD = 0.05            # Angstrom, per Cartesian component (~0.1 A displacement magnitude)
SEED = 42                    # fixed seed so every group member gets the same structure

# Cutoff scan (8-atom cell, fixed k-mesh)
ECUT_LIST = [30, 40, 50, 60, 70, 80]   # Ry
K_FOR_ECUT_SCAN = 4                    # NxNxN on the 8-atom cell

# k-point scan (64-atom production cell, fixed cutoff).
# Set ECUT_FOR_K_SCAN to the value you chose from the cutoff scan.
ECUT_FOR_K_SCAN = 40                   # Ry
K_LIST = [1, 2, 3, 4]                  # NxNxN on the 64-atom cell

# Conversions
RY_BOHR_TO_EV_ANG = 13.605693123 / 0.529177210903   # Ry/bohr -> eV/Angstrom
RY_TO_EV = 13.605693123


def build_si(n=1, rattle=True):
    """
    Builds 8-atom, Diamond-cubic Si cell with n^3 copies

    Returns: 3x3 cell in Angstrom and Cartesian positions in Angstrom
    """
    fcc = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
    basis = np.vstack([fcc, fcc + 0.25])                 # 8 atoms, fractional (cubic cell)
    frac = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                frac.append((basis + [i, j, k]) / n)
    frac = np.vstack(frac)
    cell = np.eye(3) * A_LATTICE * n
    pos = frac @ cell
    if rattle:
        rng = np.random.default_rng(SEED)
        pos = pos + rng.normal(0.0, RATTLE_STD, pos.shape)
    return cell, pos


def write_pw_input(path, prefix, cell, pos, ecutwfc, kmesh):
    """Write a pw.x single-point (scf) input that prints forces and stress."""
    nat = len(pos)
    lines = [
        "&CONTROL",
        "  calculation = 'scf'",
        f"  prefix      = '{prefix}'",
        "  outdir      = './tmp'",
        f"  pseudo_dir  = '{PSEUDO_DIR}'",
        "  tprnfor     = .true.     ! print forces",
        "  tstress     = .true.     ! print stress",
        "  disk_io     = 'none'     ! we only need the .out file",
        "/",
        "&SYSTEM",
        "  ibrav     = 0",
        f"  nat       = {nat}",
        "  ntyp      = 1",
        f"  ecutwfc   = {ecutwfc}",
        f"  ecutrho   = {ecutwfc * DUAL}",
        "  occupations = 'fixed'    ! Si is a semiconductor",
        "  nspin     = 1",
        "/",
        "&ELECTRONS",
        f"  conv_thr  = {CONV_THR}",
        "  mixing_beta = 0.5",
        "  electron_maxstep = 200",
        "/",
        "ATOMIC_SPECIES",
        f"  Si {MASS_SI} {PSEUDO_FILE}",
        "K_POINTS automatic",
        f"  {kmesh} {kmesh} {kmesh} 0 0 0",
        "CELL_PARAMETERS angstrom",
    ]
    lines += ["  %.10f %.10f %.10f" % tuple(v) for v in cell]
    lines.append("ATOMIC_POSITIONS angstrom")
    lines += ["  Si %.10f %.10f %.10f" % tuple(p) for p in pos]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def parse_pw_output(path):
    """Return dict(energy_eV, forces_eV_per_A (Nx3), stress_kbar (3x3))."""
    text = open(path).read()
    if "JOB DONE" not in text:
        raise RuntimeError(f"{path}: calculation did not finish")
    e = re.findall(r"^!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", text, re.M)
    energy = float(e[-1]) * RY_TO_EV

    forces = []
    lines = text.splitlines()
    start = next(i for i, l in enumerate(lines) if "Forces acting on atoms" in l)
    for l in lines[start + 1:]:
        m = re.match(r"\s+atom\s+\d+\s+type\s+\d+\s+force =\s+(\S+)\s+(\S+)\s+(\S+)", l)
        if m:
            forces.append([float(x) for x in m.groups()])
        elif forces:
            break
    forces = np.array(forces) * RY_BOHR_TO_EV_ANG

    s0 = next(i for i, l in enumerate(lines) if "total   stress" in l)
    stress = np.array([[float(x) for x in lines[s0 + 1 + r].split()[3:6]] for r in range(3)])
    return dict(energy=energy, forces=forces, stress=stress)
