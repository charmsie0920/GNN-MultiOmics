"""Number formatting shared by the results page, the analytics page and the charts."""

from __future__ import annotations


def format_ic50(value: float) -> str:
    """Format a predicted IC50 (uM) with precision that scales to its magnitude.

    A flat `.1f` rounds any sub-0.05 uM prediction (a real, very sensitive
    result -- GDSC ln_ic50 values span roughly -10 to 10, i.e. ~5e-5 to
    ~22000 uM) down to a misleading "0.0". Below 1 uM, show 3 significant
    figures instead so small-but-real values stay visible.
    """
    if value == 0:
        return "0.0"
    if abs(value) >= 1:
        return f"{value:.1f}"
    return f"{value:.3g}"
