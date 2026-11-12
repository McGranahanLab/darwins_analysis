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
                      (0.5 at 1, full credit at >= ALPACA_THRESHOLD). For a mutation-supported
                      GD the worst value rises to ALPACA_MUT_FLOOR (0.9) at >= MIN_SUPPORT doubled
                      mutations (Piotr, 2026-11-12, LTX1348: ALPACA, blind to PGDD, moves parallel
                      GDs in two clones onto their absent parent by parsimony; mutations win).
                      Parallel GDs under such a parent (no GD on it, absent from every sample,
                      ALPACA ratio >= ALPACA_THRESHOLD, GDs in >= 2 of its child branches) use the
                      parent's ratio instead (flag `alpaca_parent(<parent>:<ratio>)`, LTX0552-T2)
  minor_clone         a GD added by the fill's minor-clone pass, for clones left below their
                      sample's level (LTX0552-T2): 0.6 for a sliver of cells, 0.9 at >= 50% of a sample
  minority_coverage   heavy: a GD forced onto a clone that covers only part of a doubled
                      sample (relaxed fill pass), because no clone explains it without
                      over-doubling others; 0.2 at 0 coverage, rising to 0.5 at >= 0.5
                      (Piotr, 2026-11-11: place it, but penalise heavily)
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
Majority ploidy prior (Piotr, 2026-11-12, LTX0981-Tumour1; `majority_prior`):
  prior_raised        a borderline sample was raised to the level most samples clearly
                      share (clonal-GD prior); 0.6 at the edge of the band, 0.85 next to 0.5
  prior_lowered       a borderline sample was lowered to it (diploid prior, stronger);
                      0.75 at the edge of the band, 0.95 next to 0.5
  Moved samples are left out of borderline_ploidy and ambiguous_solution: their
  uncertainty is counted once, here. prior_lowered is tumour-only (no GD carries it).
All values are in `wgd.params` / `params.toml`.
GD order is deliberately not used (Piotr, 2026-11-10).

Per-event score (2026-11-11): the same factors restricted to one GD:
borderline_ploidy over the event's samples, low_mut_support of that event,
alpaca_disagrees of its clone, minority_coverage / fill_tie / multicluster_parent_child
of that event, ambiguous_solution if its clone's GD count changes in an alternative,
and reversion if a reverted clone has this clone as nearest doubled ancestor.
ploidy_path_mismatch, unplaced_events and fill_count_differs are tumour-only.
"""

from __future__ import annotations

import numpy as np

PLOIDY_MARGIN = 0.10
MIN_SUPPORT = 10
ALPACA_THRESHOLD = 1.5
ALPACA_MUT_FLOOR = 0.9
MAJORITY_COVERAGE = 0.5
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


def alpaca_factor(ratio: float, support: int | None = None) -> float:
    """alpaca_disagrees: `graded` from 1.0 to ALPACA_THRESHOLD, with its lowest value raised for a
    mutation-supported GD: 0.5 with no doubled mutations, ALPACA_MUT_FLOOR at >= MIN_SUPPORT."""
    x = float(np.clip((ratio - 1.0) / (ALPACA_THRESHOLD - 1.0), 0.0, 1.0))
    lo = 0.5
    if support is not None:
        lo = 0.5 + (ALPACA_MUT_FLOOR - 0.5) * float(np.clip(support / MIN_SUPPORT, 0.0, 1.0))
    return lo + (1 - lo) * x


def majority_prior(targets: dict[str, int], fracs: dict[str, tuple[float, float, float]],
                   band: float, min_clear_fraction: float, borderline_majority: bool = True
                   ) -> tuple[dict[str, int], dict[str, tuple[int, int, float, str]]]:
    """Majority ploidy prior, one GD level (boundary) at a time (Piotr, 2026-11-12).

    At boundary j (>= j GDs, genome fraction f[j-1] vs 0.5) a sample is *clear* when its
    fraction is at least PLOIDY_MARGIN from 0.5 and *borderline* when it is less than `band`
    from it. Moves are only between levels j-1 and j:
      clear majority (LTX0981): more than `min_clear_fraction` of all samples are clearly on
          one side; borderline samples on the other side join it.
      borderline majority (LTX0693, `borderline_majority`): no sample is clear at this
          boundary, every sample is borderline, and more than `min_clear_fraction` of them
          are on one side; the rest join it.
    Raising follows a doubled majority (clonal-GD prior); lowering follows an undoubled one
    (diploid prior when j = 1). Returns (new targets, {sample: (old, new, distance from 0.5,
    "clear" | "borderline")}).
    """
    new = dict(targets)
    moves: dict[str, tuple[int, int, float, str]] = {}
    n = len(targets)
    for j in (1, 2, 3):
        x = {s: fracs[s][j - 1] - 0.5 for s in targets
             if fracs.get(s) is not None and fracs[s][j - 1] is not None}
        if not x:
            continue
        up = [s for s, v in x.items() if v >= PLOIDY_MARGIN]
        down = [s for s, v in x.items() if v <= -PLOIDY_MARGIN]
        border_up = [s for s, v in x.items() if 0 <= v < band]
        border_down = [s for s, v in x.items() if -band < v < 0]
        if len(up) > min_clear_fraction * n:
            side, kind = "up", "clear"
        elif len(down) > min_clear_fraction * n:
            side, kind = "down", "clear"
        elif (borderline_majority and not up and not down and len(border_up) + len(border_down) == n
              and n > 0):
            if len(border_up) > min_clear_fraction * n:
                side, kind = "up", "borderline"
            elif len(border_down) > min_clear_fraction * n:
                side, kind = "down", "borderline"
            else:
                continue
        else:
            continue
        for s in (border_down if side == "up" else border_up):
            frm, to = (j - 1, j) if side == "up" else (j, j - 1)
            if new[s] == frm:
                new[s] = to
                moves[s] = (targets[s] if s not in moves else moves[s][0], to, float(abs(x[s])), kind)
    return new, moves


def prior_factor(d: float, band: float, lo_hi: tuple[float, float]) -> float:
    """`hi` for a sample right at 0.5, falling linearly to `lo` at the edge of the band."""
    lo, hi = lo_hi
    return hi - (hi - lo) * float(np.clip(d / band, 0.0, 1.0))


MINOR_CLONE_FACTOR = (0.6, 0.9)


def minor_clone(share: float) -> float:
    """A GD added only for minor clones (fill minor-clone pass, LTX0552-T2): MINOR_CLONE_FACTOR[0]
    when it brings almost no cells to target, rising to [1] at MAJORITY_COVERAGE of a sample
    (largest share over its samples)."""
    lo, hi = MINOR_CLONE_FACTOR
    return lo + (hi - lo) * float(np.clip(share / MAJORITY_COVERAGE, 0.0, 1.0))


def forced(coverage: float) -> float:
    """Heavy penalty for a forced (relaxed-pass) GD: 0.2 at 0 coverage, 0.5 at >= MAJORITY_COVERAGE."""
    return 0.2 + 0.3 * float(np.clip(coverage / MAJORITY_COVERAGE, 0.0, 1.0))


TUMOUR_ONLY = ("borderline_ploidy", "ploidy_path_mismatch", "unplaced_events", "fill_count_differs",
               "prior_lowered")


def aggregate(event_factors: list[list], tumour_factors: list) -> tuple[float, str | None]:
    """Tumour score from its events' factors plus tumour-level factors: per code, the worst
    (lowest) value, then the product. Used both for the full table and after rejecting events."""
    worst: dict[str, tuple[float, str]] = {}
    for code, f, detail in [x for fs in event_factors for x in fs] + list(tumour_factors):
        if code not in worst or f < worst[code][0]:
            worst[code] = (f, detail)
    return score([(c, f, d) for c, (f, d) in worst.items()])


def score(factors: list[tuple[str, float, str]]) -> tuple[float, str | None]:
    """(product of factors, codes for factors < 1). factors: (code, factor, detail)."""
    q = 1.0
    codes = []
    for code, f, detail in factors:
        if f < 1:
            q *= f
            codes.append(f"{code}({detail})" if detail else code)
    return round(q, 3), ";".join(codes) or None
