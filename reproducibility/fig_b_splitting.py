"""Paper B flagship figure: spin-symmetry breaking of autoionization
widths vs Zalpha, at two reduced screenings (near-threshold s=0.24 and
deep-continuum s=0.30), from f3_table.csv.

=======================================================================
SHIPPED COPY -- INPUT MAP  (feast-spectral: exploration/dirac_cs/
fig_b_splitting.py, verbatim except for the four marked lines below)
=======================================================================
This is the "ratio script": lines ~52-60 compute the width ratio
Gamma(kappa=-2)/Gamma(kappa=+1) = Gamma(2p_3/2)/Gamma(2p_1/2) that is
plotted in the LEFT panel of the flagship figure, directly from the
solver's production table.  Right panel: complex-plane trajectories at
Z = 80.

  motor file (generator)                shipped name in this deposit
  ------------------------------------  ----------------------------
  exploration/dirac_cs/f3_table.csv     ../data/solver_grid_2p.csv
  exploration/dirac_cs/f3_traj.csv      ../data/f3_traj_Z80.csv

Both shipped CSVs are byte-identical to the motor files they came from.
Columns: see ../data/README.md.  `Gamma` is the FULL width; the right
panel plots -Gamma/2, the imaginary part in the E = E_r - i*Gamma/2
convention.

Deviations from the motor original, all confined to plumbing (marked
`# SHIPPED:` inline) -- no plotting, selection or arithmetic changed:
  1. C_AU is taken from a literal instead of `from dirac import C_AU`,
     so the script runs without the solver repo checked out.
  2. The two csv.DictReader paths resolve the shipped names above,
     falling back to the motor names if run inside the solver repo.
  3. The output goes next to this script instead of <repo>/figures/.

Run:  python fig_b_splitting.py      (needs numpy + matplotlib only)
=======================================================================
"""

import csv
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# SHIPPED (1): motor line was `from dirac import C_AU`.
# c = 1/alpha in Hartree atomic units, CODATA 2022 (feast-spectral
# src/dirac.py: C_AU = 137.035999084).
C_AU = 137.035999084

HERE = Path(__file__).resolve().parent


def _input(shipped, motor):
    """SHIPPED (2): prefer the deposited CSV name, fall back to the
    motor name when this file sits in exploration/dirac_cs/."""
    for p in (HERE.parent / "data" / shipped, HERE / shipped,
              HERE / motor):
        if p.is_file():
            return p
    raise FileNotFoundError(f"{shipped} (or {motor}) not found near {HERE}")


plt.rcParams.update({
    "font.family": "serif", "font.size": 10,
    "axes.labelsize": 10.5, "xtick.labelsize": 9.5, "ytick.labelsize": 9.5,
    "legend.fontsize": 8.8, "axes.grid": True, "grid.alpha": 0.3,
    "figure.dpi": 200, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})

data = {}
for r in csv.DictReader(_input("solver_grid_2p.csv", "f3_table.csv").open()):
    if not r["Er"] or r["Er"] == "nan":
        continue
    data[(float(r["Z"]), float(r["s"]), int(r["kappa"]))] = (
        float(r["Er"]), float(r["Gamma"]))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(390.0 / 72.27, 2.9))
for s, mk, col in ((0.24, "o", "tab:red"), (0.30, "s", "tab:blue")):
    za, ratio = [], []
    for Z in (1.0, 20.0, 40.0, 60.0, 80.0, 92.0):
        k1, k2 = (Z, s, 1), (Z, s, -2)
        if k1 in data and k2 in data:
            za.append(Z / C_AU)
            ratio.append(data[k2][1] / data[k1][1])
    ax1.semilogy(za, ratio, mk + "-", color=col, ms=4.5, lw=1.1,
                 label=f"$s = {s}$")
# bound-partner region at s=0.24 (Z >= 92): ratio formally infinite
ax1.axvspan(92 / C_AU, 0.98, color="tab:red", alpha=0.12)
ax1.text(0.80, 2e2, "$2p_{1/2}$ bound,\n$2p_{3/2}$ resonant",
         fontsize=8, color="tab:red", ha="center")
ax1.set_xlabel(r"$Z\alpha$")
ax1.set_ylabel(r"$\Gamma_{2p_{3/2}}\,/\,\Gamma_{2p_{1/2}}$")
ax1.legend(loc="upper left")

# trajectories at Z=80
traj = {1: [], -2: []}
for r in csv.DictReader(_input("f3_traj_Z80.csv", "f3_traj.csv").open()):
    traj[int(r["kappa"])].append((float(r["Er"]), float(r["Gamma"])))
for kappa, name, col in ((1, r"$2p_{1/2}$", "tab:orange"),
                         (-2, r"$2p_{3/2}$", "tab:green")):
    pts = np.array(traj[kappa])
    ax2.plot(pts[:, 0], -pts[:, 1] / 2, "o-", color=col, ms=4, lw=1.1,
             label=name)
ax2.set_xlabel(r"$E_r$ (a.u.)")
ax2.set_ylabel(r"$-\Gamma/2$ (a.u.)")
ax2.set_title(r"$Z=80$, $s = 0.245 \to 0.30$", fontsize=9)
ax2.legend(loc="lower left")
fig.tight_layout(w_pad=1.4)
out = HERE / "figB_splitting"          # SHIPPED (3): motor wrote to <repo>/figures/
for ext in ("pdf", "png"):
    fig.savefig(f"{out}.{ext}")
print("figB_splitting done")

# SHIPPED: the ratio actually plotted in the left panel, printed so the
# number can be read without opening the figure.
print("\nleft panel, Gamma(2p_3/2)/Gamma(2p_1/2):")
for s in (0.24, 0.30):
    for Z in (1.0, 20.0, 40.0, 60.0, 80.0, 92.0):
        k1, k2 = (Z, s, 1), (Z, s, -2)
        if k1 in data and k2 in data:
            print(f"   s={s:.2f}  Z={Z:>5.1f}  Zalpha={Z/C_AU:.4f}"
                  f"   ratio = {data[k2][1]/data[k1][1]:.6f}")
