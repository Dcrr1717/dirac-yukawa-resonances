"""Complex-scaled radial Dirac: shared assembly (Plan B core).

The shifted complex-symmetric pencil is built from the SAME certified
real matrices as the Hermitian solver; see f0_feasibility.py for the
gate battery that validates this module.
"""

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from dirac import C_AU  # noqa: E402

C = C_AU


def cs_pencil(B, theta, lam=0.0):
    """(H_theta - c^2 S) on basis B; lam=0 pure Coulomb (exact),
    lam>0 Yukawa (screening factor at complex argument, converged
    quadrature)."""
    c, Z, r = B.c, B.Z, B.r
    ph = np.exp(-1j * theta)
    f = B.f
    ME, MEv = B._me, B._mev
    K = ME(f, 2.0 * c * B.h)
    if lam == 0.0:
        Vth = ph * (-Z) * MEv(f, 1.0 / r, f).astype(complex)
    else:
        scr = np.exp(-lam * r * np.exp(1j * theta))
        Vth = ph * (-Z) * np.einsum('ik,k,jk,k->ij', f, scr / r, f, B.W)
    n = B.nb
    Hs = np.block([[Vth, c * ph * K],
                   [c * ph * K.T, Vth - 2 * c * c * np.eye(n)]])
    return 0.5 * (Hs + Hs.T)


def resonance_eigs(Hs, re_win, im_win):
    ev = np.linalg.eigvals(Hs)
    return ev[(ev.real > re_win[0]) & (ev.real < re_win[1])
              & (ev.imag > im_win[0]) & (ev.imag < im_win[1])]
