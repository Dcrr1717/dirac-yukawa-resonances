"""Plan B, Phase F3: production tables and the spin-symmetry-breaking map.

E-B3a  Reference table for the 2p doublet (kappa = +1: 2p1/2,
       kappa = -2: 2p3/2) of the Dirac-Yukawa problem, Z = 1..130,
       reduced screening s = lambda/Z in {0.24, 0.30}.
       HONEST-DIGITS PROTOCOL:每 state evaluated on a 2x2x2 grid
       (theta pair x nb pair x beta pair); quoted error = max spread;
       Gamow-shooting anchor where Gamma <~ 4 (the double-precision
       envelope of the incoming-coefficient extraction).
E-B3b  The flagship map: fine-structure splitting of widths vs Zalpha.
       Expectation from E-A3/F2: percent-level at low Z, ORDERS OF
       MAGNITUDE at strong coupling; partner may even remain BOUND
       (s < s_c(2p1/2)) while the other resonates -- report that too.
E-B3c  Complex-plane trajectories vs s at Z = 80 for both partners.

Pre-registered: P3.1 both partners converge with spread <= 1e-8 for
narrow/moderate states; P3.2 the width ratio Gamma(3/2)/Gamma(1/2) at
fixed s grows monotonically with Zalpha.
"""

import argparse
import csv
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

# ---------------------------------------------------------------------
# OPT-IN per-point export of the honest-digits (theta x nb x beta) grid.
#
# measure() below already evaluates the full 2x2x2 grid but returns only
# the collapsed (value, spread) pair.  With --grid-dump the individual
# grid points are additionally recorded here.  When _GRID_TARGETS is
# empty -- i.e. whenever the script is run or imported WITHOUT that flag
# -- _grid_dump_on() is False, no row is ever appended, and every number
# computed, printed or written by this file is bit-for-bit what it was
# before this block existed.  Nothing in the default path is touched.
# ---------------------------------------------------------------------
_GRID_TARGETS = set()      # {(Z, kappa)} requested via --grid-dump
_GRID_ROWS = []            # collected per-(theta, nb, beta) samples


def _grid_dump_on(Z, kappa):
    return (float(Z), int(kappa)) in _GRID_TARGETS


def _grid_record(Z, kappa, lam, nb, bf, beta, thetas, pair, dev, matched,
                 val_index):
    """One (nb, beta) cell -> one row per theta actually evaluated."""
    s = lam / Z if Z else float("nan")
    for th, e in zip(thetas, pair):
        _GRID_ROWS.append(dict(
            state=f"Z{Z:g}_s{s:g}_kappa{kappa:+d}", Z=float(Z), s=s,
            kappa=int(kappa), lam=float(lam), nb=int(nb),
            beta_factor=float(bf), beta=float(beta), theta=float(th),
            Er=float(e.real), Gamma=float(-2.0 * e.imag),
            theta_pair_dev=float(dev),
            theta_matched=int(bool(matched)),
            _val_index=val_index if th == thetas[-1] else -1,
            enters_spread=0, is_reference=0, dev_from_ref=float("nan")))


def _grid_finalize(start, vals, ref, keep_mask):
    """Annotate the rows produced by one measure() call with the
    collapse bookkeeping (which grid point became the reference, which
    survived the 5% keep filter, and each one's deviation from the
    reference) so that `spread` can be recomputed from the CSV alone."""
    for row in _GRID_ROWS[start:]:
        k = row.pop("_val_index", -1)
        if k < 0 or not keep_mask[k]:
            continue
        row["enters_spread"] = 1
        row["dev_from_ref"] = float(abs(vals[k] - ref))
        row["is_reference"] = int(k == len(vals) - 1)


def measure(Z, kappa, lam, re_win=(1e-6, 200.0), im_win=(-80.0, 1e-12),
            theta_pair=(0.45, 0.60)):
    """Honest-digits resonance: 2 theta x 2 nb x 2 beta grid; returns
    (value from the largest basis, max spread, n_found).

    theta_pair is the rotation-angle pair of the theta-stability selector.
    The default (0.45, 0.60) serves narrow and moderate states; broad
    states (Gamma/E_r >~ 2, Methods sec. 'Domains of the two methods')
    require the larger pair theta=(0.65, 0.80), to which the default path
    re-selects automatically (see the tail of this function)."""
    nu = Z / C
    gam = float(np.sqrt(kappa * kappa - nu * nu))
    vals = []
    dump = _grid_dump_on(Z, kappa)               # opt-in; False by default
    gstart = len(_GRID_ROWS)
    for nb in (200, 240):
        for bf in (0.25, 0.40):
            B = RadialDiracBasis(nb, kappa, Z, bf * max(Z, 1.0), mu=gam,
                                 quad="cusp", nquad=nb + 100, cert_tol=1e-5)
            picks = []
            thetas = []
            for th in theta_pair:
                Hs = cs_pencil(B, th, lam=lam)
                win = resonance_eigs(Hs, re_win, im_win)
                if win.size:
                    picks.append(win)
                    thetas.append(th)
            if len(picks) == 2:
                d = np.abs(picks[1][None, :] - picks[0][:, None])
                i, j = np.unravel_index(np.argmin(d), d.shape)
                if d[i, j] < 1e-2 * max(abs(picks[1][j]), 1e-3):
                    vals.append(picks[1][j])
                    if dump:
                        _grid_record(Z, kappa, lam, nb, bf,
                                     bf * max(Z, 1.0), thetas,
                                     (picks[0][i], picks[1][j]),
                                     float(d[i, j]), True, len(vals) - 1)
                elif dump:
                    _grid_record(Z, kappa, lam, nb, bf, bf * max(Z, 1.0),
                                 thetas, (picks[0][i], picks[1][j]),
                                 float(d[i, j]), False, -1)
            elif dump:
                # theta pair incomplete: at least one angle found no
                # eigenvalue in the window -> the cell contributes
                # nothing; record it as a NaN placeholder so the CSV
                # still shows all eight attempted grid points.
                _grid_record(Z, kappa, lam, nb, bf, bf * max(Z, 1.0),
                             theta_pair,
                             (complex("nan+nanj"), complex("nan+nanj")),
                             float("nan"), False, -1)
    _narrow = (theta_pair == (0.45, 0.60))
    if not vals:
        if dump:
            _grid_finalize(gstart, [], None, [])
        # Broad states are under-rotated at the default pair and found at
        # neither angle; re-select once at the Methods' broad pair.
        if _narrow and not dump:
            return measure(Z, kappa, lam, re_win, im_win,
                           theta_pair=(0.65, 0.80))
        return None, np.inf, 0
    vals = np.array(vals)
    ref = vals[-1]
    keep = vals[np.abs(vals - ref) < 0.05 * max(abs(ref), 1e-3)]
    if dump:
        _grid_finalize(gstart, vals, ref,
                       np.abs(vals - ref) < 0.05 * max(abs(ref), 1e-3))
    # Broad states (Gamma/E_r >~ 2) require the larger rotation pair; at
    # the default pair they are either missed or captured as a static
    # rotated-continuum artifact.  Re-select once at theta=(0.65, 0.80).
    if (_narrow and not dump and ref.real > 0
            and -2.0 * ref.imag > 2.0 * ref.real):
        return measure(Z, kappa, lam, re_win, im_win,
                       theta_pair=(0.65, 0.80))
    return ref, float(np.max(np.abs(keep - ref))) if keep.size > 1 else np.inf, keep.size


def _grid_dump_main(targets, out_points, out_summary):
    """--grid-dump driver.  Runs measure() (and, where the width allows,
    the same Gamow-shooting anchor the production sweep uses) on the
    requested (Z, s, kappa) states ONLY, writes the per-point grid and a
    summary in the exact column layout of f3_table.csv, and exits.  The
    production tables f3_table.csv / f3_traj.csv are NEVER written in
    this mode, so an opt-in dump run cannot perturb them."""
    summary = []
    print("opt-in grid dump; production tables NOT written in this mode")
    print(f"{'Z':>5} {'s':>6} {'kap':>4} {'E_r':>18} {'Gamma':>15}"
          f" {'spread':>10} {'shoot':>10} {'n':>3}")
    for Z, s, kappa in targets:
        lam = s * Z
        tstate = time.perf_counter()
        e, spread, n = measure(Z, kappa, lam)
        sh = np.nan
        if e is not None and -2 * e.imag < 4.0 and e.real > 1e-5:
            try:                                  # same anchor as the sweep
                Rr = max(18.0 / lam, 2.0) if lam > 0 else 20.0
                es = refine(complex(e) + C**2, Z, lam, kappa, Rr) - C**2
                sh = abs(e - es)
            except Exception as exc:
                print(f"      shooting failed: {type(exc).__name__}: {exc}")
        if e is None:
            print(f"{Z:>5.0f} {s:>6.3f} {kappa:>4}   (no resonance --"
                  f" bound or dissolved)   [{time.perf_counter()-tstate:.1f} s]")
            summary.append(dict(Z=Z, s=s, kappa=kappa, Er=np.nan,
                                Gamma=np.nan, spread=np.nan, shoot=np.nan))
            continue
        print(f"{Z:>5.0f} {s:>6.3f} {kappa:>4} {e.real:>18.12f}"
              f" {-2 * e.imag:>15.9e} {spread:>10.2e} {sh:>10.2e} {n:>3}"
              f"   [{time.perf_counter()-tstate:.1f} s]")
        summary.append(dict(Z=Z, s=s, kappa=kappa, Er=e.real,
                            Gamma=-2 * e.imag, spread=spread, shoot=sh))
    if _GRID_ROWS:
        cols = [k for k in _GRID_ROWS[0] if not k.startswith("_")]
        with Path(out_points).open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            for r in _GRID_ROWS:
                w.writerow(r)
        print(f"\ngrid points -> {out_points}  ({len(_GRID_ROWS)} rows)")
    if out_summary and summary:
        with Path(out_summary).open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(summary[0].keys()))
            w.writeheader()
            for r in summary:
                w.writerow(r)
        print(f"summary     -> {out_summary}  ({len(summary)} rows)")
    print(f"\ntotal: {time.perf_counter()-t0:.1f} s")


_ap = argparse.ArgumentParser(add_help=False)
_ap.add_argument("--grid-dump", action="append", default=[], metavar="Z,s,kappa",
                 help="opt-in: export every (theta, nb, beta) grid point for "
                      "this state; repeatable. Implies dump-only mode.")
_ap.add_argument("--grid-dump-out", default=None, metavar="PATH")
_ap.add_argument("--grid-summary-out", default=None, metavar="PATH")
_args, _rest = _ap.parse_known_args()
if _args.grid_dump:
    _targets = []
    for _spec in _args.grid_dump:
        _Z, _s, _k = _spec.split(",")
        _targets.append((float(_Z), float(_s), int(_k)))
        _GRID_TARGETS.add((float(_Z), int(_k)))
    _grid_dump_main(
        _targets,
        _args.grid_dump_out or Path(__file__).with_name("f3_grid_dump.csv"),
        _args.grid_summary_out)
    sys.exit(0)
# ---- default path resumes here; unchanged from here down -------------


rows = []
print("E-B3a/b production (2p doublet), honest-digits protocol:")
print(f"{'Z':>4} {'s':>5} {'kap':>4} {'E_r':>16} {'Gamma':>13}"
      f" {'spread':>9} {'shoot':>9}")
for Z in (1.0, 20.0, 40.0, 60.0, 80.0, 92.0, 110.0, 130.0):
    nu = Z / C
    for s in (0.24, 0.30):
        lam = s * Z
        pair = {}
        for kappa in (1, -2):
            e, spread, n = measure(Z, kappa, lam)
            tag = ""
            sh = np.nan
            if e is not None and -2 * e.imag < 4.0 and e.real > 1e-5:
                try:
                    Rr = max(18.0 / lam, 2.0) if lam > 0 else 20.0
                    es = refine(complex(e) + C**2, Z, lam, kappa, Rr) - C**2
                    sh = abs(e - es)
                except Exception:
                    pass
            if e is None:
                print(f"{Z:>4.0f} {s:>5.2f} {kappa:>4}   (no resonance --"
                      f" bound or dissolved)")
                rows.append(dict(Z=Z, s=s, kappa=kappa, Er=np.nan,
                                 Gamma=np.nan, spread=np.nan, shoot=np.nan))
                continue
            pair[kappa] = e
            print(f"{Z:>4.0f} {s:>5.2f} {kappa:>4} {e.real:>16.9f}"
                  f" {-2 * e.imag:>13.6e} {spread:>9.1e} {sh:>9.1e}")
            rows.append(dict(Z=Z, s=s, kappa=kappa, Er=e.real,
                             Gamma=-2 * e.imag, spread=spread, shoot=sh))
        if len(pair) == 2:
            g1, g2 = -2 * pair[1].imag, -2 * pair[-2].imag
            print(f"{'':>14} doublet: dE_r = {pair[-2].real - pair[1].real:+.4e}"
                  f"   Gamma(3/2)/Gamma(1/2) = {g2 / g1:.3f}")

out = Path(__file__).with_name("f3_table.csv")
with out.open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    for r in rows:
        w.writerow(r)

# ---------------- E-B3c: trajectories at Z=80 ---------------------------
print("\nE-B3c trajectories (Z=80, s sweep):")
traj = {1: [], -2: []}
for s in (0.245, 0.26, 0.28, 0.30, 0.33, 0.36):
    for kappa in (1, -2):
        e, spread, _ = measure(80.0, kappa, s * 80.0)
        if e is not None:
            traj[kappa].append((s, e.real, -2 * e.imag))
for kappa, name in ((1, "2p1/2"), (-2, "2p3/2")):
    print(f"   {name}: " + "  ".join(
        f"s={s:.3f}:({er:.3f},{ga:.3f})" for s, er, ga in traj[kappa]))

with Path(__file__).with_name("f3_traj.csv").open("w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["kappa", "s", "Er", "Gamma"])
    for kappa, pts in traj.items():
        for s, er, ga in pts:
            w.writerow([kappa, s, er, ga])
print(f"\ntotal: {time.perf_counter()-t0:.1f} s")
