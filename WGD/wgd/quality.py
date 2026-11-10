"""Per-tumour quality score in [0, 1] (Piotr, 2026-11-10).

The score is the product of factors in (0, 1], one per piece of evidence. Every
factor below 1 adds a code to `quality_code` (`;`-joined, with the value that
triggered it). 1 = no concern. Threshold it later to get High/Low back.

Graded factors run from 0.5 (worst) to 1 (at or past the threshold):
  borderline_ploidy   smallest distance of any sample's deciding genome fraction from 0.5
                      (full credit at >= PLOIDY_MARGIN)
  low_mut_support     smallest supporting doubled-mutation count of a mutation-supported GD
                      (full credit at >= MIN_SUPPORT)
  alpaca_disagrees    smallest ALPACA max(ratio_A, ratio_B) among clones carrying a GD
                      (0.5 at 1, full credit at >= ALPACA_THRESHOLD)
Fixed factors:
  ambiguous_solution  0.5   the placement changes when borderline samples are flipped
                            or the clone-majority cut-off moves to 0.4/0.6
  reversion           0.2   a clone below a doubled clone looks diploid again
                            (median ALPACA total CN <= 2)
  ploidy_path_mismatch 0.5  GDs on a root-leaf path don't match the samples' num_gds (H)
  unplaced_events     0.5   some PGDD events could not be put on an ALPACA clone
  multicluster_parent_child 0.7
  fill_tie            0.8   a fill GD was chosen by the ALPACA score alone
  fill_count_differs  0.9   the fill added a different number of GDs than PGDD's ploidy-only events
GD order is deliberately not used (Piotr, 2026-11-10).
"""

from __future__ import annotations

import numpy as np

PLOIDY_MARGIN = 0.10
MIN_SUPPORT = 10
ALPACA_THRESHOLD = 1.5
FIXED = {
    "ambiguous_solution": 0.5,
    "reversion": 0.2,
    "ploidy_path_mismatch": 0.5,
    "unplaced_events": 0.5,
    "multicluster_parent_child": 0.7,
    "fill_tie": 0.8,
    "fill_count_differs": 0.9,
}


def graded(value: float, full_at: float, zero_at: float = 0.0) -> float:
    """0.5 at `zero_at`, rising linearly to 1 at `full_at`."""
    x = (value - zero_at) / (full_at - zero_at)
    return 0.5 + 0.5 * float(np.clip(x, 0.0, 1.0))


def ploidy_margins(targets: dict[str, int], fracs: dict[str, tuple[float, float, float]]) -> dict[str, float]:
    """Per sample: distance of the deciding genome fraction(s) from 0.5.

    `fracs[s]` = (genome_frac_nMaj2, nMaj3, nMaj4), i.e. the fraction supporting >= 1, 2, 3 GDs.
    A sample with n GDs is decided by fraction n (must be > 0.5) and fraction n+1 (must not be).
    """
    out = {}
    for s, n in targets.items():
        f = fracs.get(s)
        if f is None:
            continue
        d = []
        if 1 <= n <= 3 and f[n - 1] is not None:
            d.append(f[n - 1] - 0.5)
        if n <= 2 and f[n] is not None:
            d.append(0.5 - f[n])
        if d:
            out[s] = min(d)
    return out


def flips(targets: dict[str, int], fracs: dict[str, tuple[float, float, float]],
          margins: dict[str, float]) -> list[dict[str, int]]:
    """Alternative targets: each borderline sample moved across its nearest boundary, then all at once."""
    moves = {}
    for s, m in margins.items():
        if m >= PLOIDY_MARGIN:
            continue
        n, f = targets[s], fracs[s]
        down = f[n - 1] - 0.5 if n >= 1 else np.inf   # how far above the boundary that gives it n GDs
        up = 0.5 - f[n] if n <= 2 else np.inf          # how far below the boundary for n + 1
        moves[s] = n - 1 if down <= up else n + 1
    alts = [{**targets, s: v} for s, v in moves.items()]
    if len(moves) > 1:
        alts.append({**targets, **moves})
    return alts


def score(factors: list[tuple[str, float, str]]) -> tuple[float, str | None]:
    """(product of factors, codes for factors < 1). factors: (code, factor, detail)."""
    q = 1.0
    codes = []
    for code, f, detail in factors:
        if f < 1:
            q *= f
            codes.append(f"{code}({detail})" if detail else code)
    return round(q, 3), ";".join(codes) or None
