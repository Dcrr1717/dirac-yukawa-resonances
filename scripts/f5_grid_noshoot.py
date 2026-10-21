"""Plan B, Phase F5: memory-lean driver for f3_production's --grid-dump.

Identical to running

    python f3_production.py --grid-dump ... --grid-dump-out ... --grid-summary-out ...

except that the module-level self-test of f2_shooting (twenty extra
basis solves at nb = 200/240, run merely as a side effect of the
`from f2_shooting import refine` line) is skipped by pre-registering a
stub module.  Nothing in the spectral path changes: measure() and
_grid_dump_main() are the untouched originals, so the eigenvalues,
spreads and CSV layout are bit-for-bit those of f3_production.

Consequence: the inline Gamow anchor is unavailable in this mode and the
`shoot` column comes out nan.  The independent shooting cross-check is
run separately, and more thoroughly (with matching-radius stability), by
f4_d_shooting_check.py.

Usage: same flags as f3_production.py --grid-dump.
"""

import sys
import types
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

_stub = types.ModuleType("f2_shooting")
_stub.refine = None          # any call raises -> f3 prints "shooting failed"
sys.modules["f2_shooting"] = _stub

if not any(a.startswith("--grid-dump") for a in sys.argv[1:]):
    sys.exit("f5_grid_noshoot.py: pass the same --grid-dump flags as f3_production.py")

import f3_production  # noqa: E402,F401  (does the work at import, then exits)
