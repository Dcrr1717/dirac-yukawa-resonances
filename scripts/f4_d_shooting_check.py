"""Plan B, Phase F4: independent Gamow-shooting cross-check of the 3d
doublet resonances, with matching-radius (R) stability.

Mirrors exactly what f2_shooting.py's CASES block does for the 2p
doublet: take the theta-stable spectral value as the secant seed, drive
the incoming-wave (Siegert) coefficient to zero on the radial Dirac
system matched to the exact free outgoing spinor, and repeat at several
matching radii R so the R-drift can be quoted alongside the deviation.

Read-only with respect to every existing artifact: it consumes the
summary written by

    python f3_production.py --grid-dump ... --grid-summary-out d_doublet_summary.csv

and writes only its own output file.
"""

import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from f2_shooting import C, refine  # noqa: E402  (module self-test runs on import)

ap = argparse.ArgumentParser()
ap.add_argument("--summary", action="append", default=[],
                help="d_doublet_summary*.csv produced by f3_production.py")
ap.add_argument("--out", default="d_doublet_shooting.csv")
ap.add_argument("--gamma-max", type=float, default=6.0,
                help="skip states broader than this (shooting domain)")
args = ap.parse_args()

rows = []
for path in args.summary or ["d_doublet_summary.csv"]:
    with Path(path).open() as fh:
        rows.extend(list(csv.DictReader(fh)))

t0 = time.perf_counter()
print("F4: 3d doublet, spectral vs INDEPENDENT Gamow shooting")
print(f"{'Z':>5} {'s':>7} {'kap':>4} {'Gamma':>12} {'R':>7} "
      f"{'rel.dev':>10} {'R-drift':>10}")

out = []
for r in rows:
    Z, s, kappa = float(r["Z"]), float(r["s"]), int(r["kappa"])
    Er, Ga = float(r["Er"]), float(r["Gamma"])
    if not np.isfinite(Er) or Ga > args.gamma_max:
        print(f"{Z:>5.0f} {s:>7.4f} {kappa:>4} {Ga:>12.4e}"
              f"       ---        ---        ---   (outside shooting domain)")
        out.append(dict(Z=Z, s=s, kappa=kappa, Er_spectral=Er, Gamma=Ga,
                        R=np.nan, Er_shoot=np.nan, Gamma_shoot=np.nan,
                        rel_dev=np.nan, R_drift=np.nan))
        continue
    lam = s * Z
    e = complex(Er, -0.5 * Ga)
    Rbase = max(18.0 / lam, 2.0) if lam > 0 else 20.0
    vals = {}
    for fac in (0.8, 1.0, 1.25, 1.5):
        R = fac * Rbase
        try:
            vals[R] = refine(e + C**2, Z, lam, kappa, R) - C**2
        except Exception as exc:                       # noqa: BLE001
            print(f"      R={R:.2f} failed: {type(exc).__name__}: {exc}")
    if not vals:
        continue
    es = vals[Rbase] if Rbase in vals else list(vals.values())[0]
    drift = max(abs(v - es) for v in vals.values()) if len(vals) > 1 else np.nan
    dev = abs(e - es) / abs(es)
    print(f"{Z:>5.0f} {s:>7.4f} {kappa:>4} {Ga:>12.4e} {Rbase:>7.2f} "
          f"{dev:>10.2e} {drift:>10.2e}")
    print(f"        spectral {e.real:.12f}{e.imag:+.12f}i")
    print(f"        shooting {es.real:.12f}{es.imag:+.12f}i")
    out.append(dict(Z=Z, s=s, kappa=kappa, Er_spectral=e.real, Gamma=Ga,
                    R=Rbase, Er_shoot=es.real, Gamma_shoot=-2 * es.imag,
                    rel_dev=dev, R_drift=drift))

with Path(args.out).open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
    w.writeheader()
    for r in out:
        w.writerow(r)
print(f"\n-> {args.out}  ({len(out)} rows)")
print(f"total: {time.perf_counter() - t0:.1f} s")
