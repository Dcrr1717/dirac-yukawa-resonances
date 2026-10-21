"""Plan B, Phase F2: INDEPENDENT validation by complex outgoing-wave
(Siegert/Gamow) shooting on the radial Dirac ODE, plus the lambda->
lambda_c consistency check against the E-A3 map.

Method (fully independent of the basis solver):
  F' = -(kappa/r) F + [(E - V + c^2)/c] G
  G' = +(kappa/r) G - [(E - V - c^2)/c] F,   V = -Z e^{-lam r}/r.
Regular start at r0 = 1e-8 with the exact indicial direction
(F,G) ~ r^gamma (1, (gamma+kappa)/nu) plus the first-order series
correction (relative start error O((2c r0)^2) ~ 1e-11); integration in
s = ln r up to r = 1 (removes the 1/r stiffness), then in r out to R.
Outward integration at complex E favors the OUTGOING (growing) Gamow
component, so the resonance condition is the vanishing of the INCOMING
coefficient, read off as the Wronskian of the regular solution with the
free outgoing spinor  out(r) = e^{ipr} (1, ipc/(E+c^2)),
p = sqrt(E^2-c^4)/c (Re p > 0):
  W(E) = F_reg G_out - G_reg F_out  -> 0.
Complex secant iteration from the basis seed.

PRE-REGISTERED PREDICTIONS:
P-F2.1 basis and shooting agree to >= 10 digits (gate G2) on >= 3 cases.
P-F2.2 W(E) is R-stable: moving R by ~20% changes E by < 1e-11.
P-F2.3 lambda scan across lambda_c (Z=40, kappa=-2): the resonance is
       born at the E-A3 threshold; E_r ~ (lam - lam_c) (ell=1 law) and
       Gamma follows the Wigner p-wave law Gamma ~ E_r^{3/2}.
"""

import sys
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "exploration" / "dirac_cs"))

from dirac import C_AU, RadialDiracBasis  # noqa: E402
from cs_core import cs_pencil, resonance_eigs  # noqa: E402

C = C_AU
t0 = time.perf_counter()


def shoot_W(E, Z, lam, kappa, R, r0=1e-8, rtol=1e-12):
    """Wronskian of the regular solution with the outgoing spinor at R."""
    nu = Z / C
    gam = np.sqrt(kappa * kappa - nu * nu)

    def V(r):
        return -Z * np.exp(-lam * r) / r

    # indicial start + first-order series (Yukawa: V = -Z/r + Z lam + ...)
    a0, b0 = 1.0, (gam + kappa) / nu
    V1 = Z * lam
    A = np.array([[gam + 1 + kappa, -nu],
                  [nu, gam + 1 - kappa]], dtype=complex)
    rhs = np.array([(E + C**2 - V1) / C * b0,
                    -(E - C**2 - V1) / C * a0])
    a1, b1 = np.linalg.solve(A, rhs)
    y0 = np.array([a0 + a1 * r0, b0 + b1 * r0], dtype=complex)

    def rhs_s(s, y):                       # s = ln r (kills 1/r stiffness)
        r = np.exp(s)
        F, G = y
        dF = -kappa * F + r * (E - V(r) + C**2) / C * G
        dG = +kappa * G - r * (E - V(r) - C**2) / C * F
        return [dF, dG]

    sol = solve_ivp(rhs_s, [np.log(r0), 0.0], y0, method='DOP853',
                    rtol=rtol, atol=1e-30)
    y1 = sol.y[:, -1]
    y1 = y1 / np.max(np.abs(y1))

    def rhs_r(r, y):
        F, G = y
        dF = -kappa / r * F + (E - V(r) + C**2) / C * G
        dG = +kappa / r * G - (E - V(r) - C**2) / C * F
        return [dF, dG]

    sol = solve_ivp(rhs_r, [1.0, R], y1, method='DOP853',
                    rtol=rtol, atol=1e-30)
    F, G = sol.y[:, -1]
    p = np.sqrt(E * E - C**4) / C
    if p.real < 0:
        p = -p
    # exact free outgoing spinor INCLUDING the centrifugal term:
    # upper F_out = Riccati-Hankel hhat_l(pr) (l of the upper component),
    # lower from the ODE row itself: G = c (F' + kappa F / r)/(E + c^2).
    ell = kappa if kappa > 0 else -kappa - 1
    z = p * R
    hh = [-1j * np.exp(1j * z), -(1 + 1j / z) * np.exp(1j * z)]
    for l in range(1, ell + 1):
        hh.append((2 * l + 1) / z * hh[l] - hh[l - 1])
    Fo = hh[ell]
    dFo = p * (hh[ell - 1] - ell / z * hh[ell]) if ell >= 1 \
        else p * np.exp(1j * z)                      # d/dr hhat_0(pr)
    Go = C * (dFo + kappa / R * Fo) / (E + C**2)
    W = F * Go - G * Fo
    return W / max(np.max(np.abs([F, G])), np.max(np.abs([Fo, Go])))


def refine(E0, Z, lam, kappa, R, tol=1e-13, itmax=40):
    """Complex secant iteration on W(E)."""
    E1 = E0 * (1 + 1e-7)
    W0, W1 = shoot_W(E0, Z, lam, kappa, R), shoot_W(E1, Z, lam, kappa, R)
    for _ in range(itmax):
        E2 = E1 - W1 * (E1 - E0) / (W1 - W0)
        if abs(E2 - E1) < tol * max(1.0, abs(E2)):
            return E2
        E0, W0, E1, W1 = E1, W1, E2, shoot_W(E2, Z, lam, kappa, R)
    return E1


def theta_stable(Z, kappa, lam, re_win, im_win, nb=200, beta=None,
                 thetas=(0.30, 0.40)):
    """Resonance selector by THETA-STABILITY: eigenvalues of the rotated
    continuum move with theta; the resonance does not.  Returns the
    window eigenvalue with the smallest displacement between the two
    rotation angles (and that displacement as a quality tag)."""
    nu = Z / C
    gam = float(np.sqrt(kappa * kappa - nu * nu))
    beta = beta or 0.25 * max(Z, 1.0)
    B = RadialDiracBasis(nb, kappa, Z, beta, mu=gam, quad="cusp",
                         nquad=nb + 120, cert_tol=1e-6)
    wins = []
    for th in thetas:
        Hs = cs_pencil(B, th, lam=lam)
        wins.append(resonance_eigs(Hs, re_win, im_win))
    if not wins[0].size or not wins[1].size:
        return None, np.inf
    d = np.abs(wins[1][None, :] - wins[0][:, None])
    i, j = np.unravel_index(np.argmin(d), d.shape)
    return wins[1][j], d[i, j]


CASES = [
    # (Z, lam, kappa, window, R)
    (1.0, 0.25, 1, ((1e-3, 0.02), (-0.05, -1e-5)), 120.0),
    (1.0, 0.25, -2, ((1e-3, 0.02), (-0.05, -1e-5)), 120.0),
    (40.0, 9.6, -2, ((0.05, 40.0), (-20.0, -1e-6)), 5.0),
    (40.0, 9.6, 1, ((0.05, 40.0), (-20.0, -1e-6)), 5.0),
    (80.0, 19.2, -2, ((0.1, 80.0), (-40.0, -1e-6)), 3.0),
    (80.0, 19.2, 1, ((0.1, 80.0), (-40.0, -1e-6)), 3.0),
]

print("P-F2.1 basis (theta-stable pick) vs independent shooting:")
for Z, lam, kappa, (re_w, im_w), R in CASES:
    Eb, qual = theta_stable(Z, kappa, lam, re_w, im_w)
    if Eb is None:
        print(f"   Z={Z:.0f} kap={kappa:+d}: no theta-stable eigenvalue")
        continue
    Es = refine(complex(Eb) + C**2, Z, lam, kappa, R) - C**2
    dev = abs(Eb - Es) / abs(Es)
    Es2 = refine(complex(Eb) + C**2, Z, lam, kappa, 1.25 * R) - C**2
    print(f"   Z={Z:>3.0f} kap={kappa:+d}: basis {Eb:.12f} (th-stab {qual:.0e})")
    print(f"                shoot {Es:.12f}   |dev| = {dev:.1e}"
          f"   R-drift = {abs(Es - Es2):.1e}")

# ---------------- P-F2.3: birth of the resonance at lambda_c ------------
print("\nP-F2.3 birth at threshold (Z=40, kappa=-2; E-A3: s_c ~ 0.2160):")
Z, kappa = 40.0, -2
rows = []
for s in (0.221, 0.225, 0.230, 0.236):
    e, qual = theta_stable(Z, kappa, s * Z, (1e-4, 3.0), (-2.0, -1e-10),
                           nb=240, beta=6.0, thetas=(0.25, 0.35))
    if e is not None and qual < 1e-3 * max(abs(e), 0.01):
        rows.append((s, e.real, -2 * e.imag))
        print(f"   s={s:.3f}: E_r = {e.real:.6f}   Gamma = {-2*e.imag:.3e}"
              f"   (th-stab {qual:.0e})")
    else:
        print(f"   s={s:.3f}: no clean theta-stable resonance"
              + (f" (best qual {qual:.0e})" if e is not None else ""))
if len(rows) >= 3:
    s_, er_, ga_ = np.array(rows).T
    lc_fit = np.polyfit(er_, s_, 1)[1]           # E_r -> 0 intercept
    wig = np.polyfit(np.log(er_), np.log(ga_), 1)[0]
    print(f"   E_r->0 extrapolated s_c = {lc_fit:.4f} (E-A3: ~0.2160)")
    print(f"   Wigner exponent d ln Gamma / d ln E_r = {wig:.3f}"
          f"  (p-wave law: 1.5)")
print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
