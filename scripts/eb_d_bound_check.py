"""Plan B: is the missing 3d partner BOUND or merely unresolved?

At Z = 130 the complex-scaled sweep finds a clean resonance for
kappa = -3 (3d_5/2) at s = 0.095 but NO theta-stable eigenvalue for
kappa = +2 (3d_3/2).  The manuscript's operational criterion for calling
a partner "bound" (rather than "we failed to find it") is the Hermitian
real-axis solver of the same certified discretization returning a
negative binding energy at that (Z, s, kappa).  This script applies
exactly that criterion, by calling ea3_critical_screening.ground_binding
-- the same function that produced the critical screenings -- on the
points in question, over the same (nb, beta_factor) grid the E-A3 sweep
uses, so a single badly scaled basis cannot fake either verdict.

Additive: writes nothing except its own scratch import CSV.
"""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# ea3 does its sweep at module level.  Give it a deliberately small but
# VALID opt-in run (one channel, one Z, a scan that really brackets that
# channel) so the import is cheap and still succeeds; we only want the
# ground_binding function it defines.
sys.argv = [sys.argv[0], "--kappas", "2", "--zs", "130",
            "--s-scan", "0.070,0.099,0.002",
            "--out", str(HERE / "_scratch_import_ea3.csv")]
import ea3_critical_screening as ea3  # noqa: E402

print("\n" + "=" * 68)
print("BOUND-vs-UNRESOLVED verdict, Hermitian real-axis solver")
print("criterion: E_b < 0  =>  bound;  E_b >= 0 / nan  =>  not bound")
print("=" * 68)
print(f"{'Z':>5} {'kap':>4} {'s':>8} {'nb':>4} {'beta':>7} {'E_b':>18}")

CASES = [(130.0, 2, 0.0925), (130.0, 2, 0.095), (130.0, 2, 0.10),
         (130.0, -3, 0.0925), (130.0, -3, 0.095), (130.0, -3, 0.10)]

verdict = {}
for Z, kappa, s in CASES:
    lam = s * Z
    ebs = []
    for nb in (90, 120):
        for bf in (0.20, 0.32):
            beta = bf * Z
            try:
                Eb = ea3.ground_binding(Z, kappa, lam, nb, beta)
            except Exception as exc:                       # noqa: BLE001
                print(f"{Z:>5.0f} {kappa:>4} {s:>8.4f} {nb:>4} {beta:>7.1f}"
                      f"   FAILED {type(exc).__name__}: {exc}")
                continue
            ebs.append(Eb)
            print(f"{Z:>5.0f} {kappa:>4} {s:>8.4f} {nb:>4} {beta:>7.1f}"
                  f" {Eb:>18.9f}")
    fin = [e for e in ebs if np.isfinite(e)]
    bound = bool(fin) and all(e < -1e-9 for e in fin)
    verdict[(Z, kappa, s)] = (bound, fin)
    tag = "BOUND" if bound else ("not bound" if fin else "no state returned")
    print(f"{'':>5} {'':>4} {'':>8}  --> {tag}"
          f"   (E_b range {min(fin):.6f} .. {max(fin):.6f})" if fin
          else f"{'':>5} {'':>4} {'':>8}  --> {tag}")
    print()

print("=" * 68)
for (Z, kappa, s), (bound, fin) in verdict.items():
    print(f"Z={Z:.0f} kappa={kappa:+d} s={s:.4f}: "
          f"{'BOUND' if bound else 'NOT bound'}")

# ---------------------------------------------------------------------
# Direct bracket of s_c at Z = 130.  The E_b -> 0 extrapolation used by
# ea3 carries a residual low bias (visible as the coarse-vs-fine scan
# shift), so bound/unbound is settled here by the solver itself: the
# last s with E_b < 0 and the first s with no state bracket s_c.
# ---------------------------------------------------------------------
print("\n" + "=" * 68)
print("DIRECT s_c BRACKET AT Z = 130 (Hermitian solver, nb=120, beta=26)")
print("=" * 68)
for kappa, grid in ((2, [0.0940, 0.0950, 0.0955, 0.0960, 0.0965, 0.0970, 0.0980]),
                    (-3, [0.0925, 0.0930, 0.0935, 0.0940, 0.0945, 0.0950])):
    lastb, firstu = None, None
    print(f"  kappa = {kappa:+d}")
    for s in grid:
        Eb = ea3.ground_binding(130.0, kappa, s * 130.0, 120, 26.0)
        ok = np.isfinite(Eb) and Eb < -1e-9
        print(f"     s={s:.4f}  E_b = {Eb:>16.9f}   {'bound' if ok else 'unbound'}")
        if ok:
            lastb = s
        elif firstu is None:
            firstu = s
    print(f"     -> s_c(Z=130, kappa={kappa:+d}) in ({lastb}, {firstu})\n")
