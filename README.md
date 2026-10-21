# Paper B — Dirac–Yukawa resonances: double-method benchmarks and threshold spin-symmetry breaking

Reproducibility repository for the Paper B manuscript (C. D. Reinoso, 2026),
prepared for The European Physical Journal D.

## Contents
- `main.tex`, `refs.bib` — manuscript source (LaTeX, `svjour` class; the
  required `svjour.cls` and `svepj.clo` are included).
- `figB_splitting.pdf` — figure.
- `data/` — tabulated results (double-method table, Shi adjudication, 2p and 3d
  doublet grids and thresholds, θ×N×β stability grid).
- `reproducibility/` — figure regenerator (`fig_b_splitting.py`), the Shi
  reproduction output, and format notes.
- `scripts/` — the production scripts: complex-scaled assembly and θ-stability
  selector (`cs_core.py`), independent Gamow shooting (`f2_shooting.py`),
  production tables (`f3_production.py`, `f5_grid_noshoot.py`), shooting
  cross-check (`f4_d_shooting_check.py`), Shi reproduction
  (`eb2c_shi_reproduction.py`), width-ratio postprocessing
  (`fig_b_splitting.py`), and the 3d critical-screening / bound-partner helpers
  (`ea3_critical_screening.py`, `eb_d_bound_check.py`).

## How to reproduce
1. Create a Python environment (NumPy/SciPy/mpmath/matplotlib).
2. Install the spectral engine (see below) so that `scripts/` can import it.
3. Run the scripts in `scripts/` to regenerate the `data/` tables and the
   figure; rebuild the PDF with `latexmk -pdf main.tex`.

## Spectral engine
The scripts in `scripts/` depend on the **spectral engine**, which is
maintained as a **separate repository** and must be **cited by its archived
DOI**. Install it into your environment (or place it on `PYTHONPATH`) before
running the scripts; it is deliberately **not** vendored into this repository.
