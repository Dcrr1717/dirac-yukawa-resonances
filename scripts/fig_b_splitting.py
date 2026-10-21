"""Paper B flagship figure: spin-symmetry breaking of autoionization
widths vs Zalpha, at two reduced screenings (near-threshold s=0.24 and
deep-continuum s=0.30), from f3_table.csv."""

import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from dirac import C_AU  # noqa: E402

plt.rcParams.update({
    "font.family": "serif", "font.size": 10,
    "axes.labelsize": 10.5, "xtick.labelsize": 9.5, "ytick.labelsize": 9.5,
    "legend.fontsize": 8.8, "axes.grid": True, "grid.alpha": 0.3,
    "figure.dpi": 200, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})

data = {}
for r in csv.DictReader((Path(__file__).parent / "f3_table.csv").open()):
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
for r in csv.DictReader((Path(__file__).parent / "f3_traj.csv").open()):
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
out = ROOT / "paper_teorema" / ".." / "figures" / "figB_splitting"
out = ROOT / "figures" / "figB_splitting"
for ext in ("pdf", "png"):
    fig.savefig(f"{out}.{ext}")
print("figB_splitting done")
