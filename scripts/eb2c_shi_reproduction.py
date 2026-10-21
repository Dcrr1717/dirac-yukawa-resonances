"""E-B2c: reproduce (and extend) the Dirac-Yukawa resonance tables of
Shi, Liu, Niu & Guo, Few-Body Syst. 55, 135 (2014).

Their setup (extracted from the paper): V = -A e^{-alpha r}/r,
hbar = M = 1, c = physical value (their theta = 1.0 rad, HO basis).
Tables: spin doublets at high l -- (kappa=+5, -6) [l=5] and
(kappa=+10, -11) [l=10] -- with strong couplings A = 110..200
(A/c = 0.8..1.46).  They report 4 digits; we bring the certified
double-method machinery (theta-stable gamma-basis + Gamow shooting).

PRE-REGISTERED PREDICTIONS:
P1 Our values agree with their 4 printed digits within their rounding
   (|dev| <~ 1e-3 absolute) for the narrow-to-moderate resonances.
P2 Basis and shooting agree to >= 9 digits on every reproduced state
   (the G2 standard), turning their 4-digit entries into benchmarks.
"""

import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "exploration" / "dirac_cs"))

from dirac import C_AU, RadialDiracBasis  # noqa: E402
from cs_core import cs_pencil, resonance_eigs  # noqa: E402
from f2_shooting import refine  # noqa: E402

C = C_AU
t0 = time.perf_counter()

# (kappa, alpha, A, [Shi's relativistic E_r - i Ei entries])
SHI = [
    (5, 2.0, 110.0, [2.2834 - 0.0549j, 3.0630 - 3.4908j]),
    (-6, 2.0, 110.0, [2.3398 - 0.0631j, 3.0454 - 3.5368j]),
    (5, 2.0, 170.0, [0.4648 - 0.0066j, 2.3913 - 2.0708j]),
    (-6, 2.0, 170.0, [0.6613 - 0.0107j, 2.3385 - 2.1687j]),
    (10, 1.0, 200.0, [0.8906 - 0.0000j, 3.2473 - 0.5827j]),
    (-11, 1.0, 200.0, [0.9312 - 0.0000j, 3.2517 - 0.5956j]),
]


def find_state(kappa, alpha, A, target, nb=240, beta=2.0):
    """Theta-stable eigenvalue nearest the (4-digit) target."""
    nu = A / C
    gam = float(np.sqrt(kappa * kappa - nu * nu))
    B = RadialDiracBasis(nb, kappa, A, beta, mu=gam, quad="cusp",
                         nquad=nb + 100, cert_tol=1e-5)
    # rotation must exceed half the resonance argument; add margin
    arg = abs(np.angle(target)) if abs(target) else 0.1
    th1 = min(max(0.30, 0.65 * arg), 1.0)
    th2 = min(th1 + 0.12, 1.15)
    picks = []
    for th in (th1, th2):
        Hs = cs_pencil(B, th, lam=alpha)
        win = resonance_eigs(Hs, (0.05, 8.0), (-8.0, 1e-12))
        if win.size:
            picks.append(win[np.argmin(np.abs(win - target))])
    if len(picks) < 2:
        return None, np.inf
    return picks[1], abs(picks[1] - picks[0])


print("Reproduction of Shi et al. (2014), Tables 1-2"
      " (their digits: 4; ours: basis + shooting):")
print(f"{'kap':>4} {'A':>4} {'state (Shi)':>19} {'ours (basis)':>26}"
      f" {'|dev|':>8} {'shoot agree':>11}")
for kappa, alpha, A, states in SHI:
    for tgt in states:
        e, q = find_state(kappa, alpha, A, tgt)
        if e is None:
            print(f"{kappa:>4} {A:>4.0f} {tgt} : not found")
            continue
        dev = abs(e - tgt)
        # independent shooting refinement (R from alpha tail)
        Rr = 18.0 / alpha
        try:
            es = refine(complex(e) + C**2, A, alpha, kappa, Rr) - C**2
            agree = abs(e - es) / abs(es)
        except Exception:
            es, agree = np.nan, np.nan
        print(f"{kappa:>4} {A:>4.0f} {tgt.real:8.4f}{tgt.imag:+8.4f}j"
              f" {e.real:13.9f}{e.imag:+13.9f}j {dev:8.1e} {agree:11.1e}")

print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
