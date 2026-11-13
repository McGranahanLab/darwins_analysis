"""Tunable parameters for step 04 (placement and quality), in one place (Piotr, 2026-11-12).

Defaults live here; `WGD/params.toml` mirrors them and is what step 04 reads
(`--params`, default `params.toml` next to `pixi.toml`). A key missing from the file keeps its
default; an unknown key is an error. `apply()` sets the module-level values in
`wgd.placement` and `wgd.quality`, so every function reads the run's values.
`place_events` applies the defaults when it is given no `Params`, so runs don't leak into
each other.
"""

from __future__ import annotations

import tomllib
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path


@dataclass(frozen=True)
class Params:
    # ploidy fill (wgd.placement)
    majority: float = 0.5                  # a clone's subtree carries a sample when > this fraction of it
    min_coverage: float = 0.1              # relaxed fill pass: smallest share of a sample a clone must cover
    fill_minor_clones: bool = True         # minor-clone pass: also double clones left below their sample's level
    minor_min_cells: float = 0.0           # minor-clone pass: a GD must bring more than this share of cells to target
    alpaca_wgd_threshold: float = 1.5      # notebook rule: alpaca_called tumour resolved iff every ratio > this

    # quality (wgd.quality)
    ploidy_margin: float = 0.10            # borderline_ploidy: full credit at this distance from 0.5
    min_support: int = 10                  # low_mut_support: full credit at this many doubled mutations
    alpaca_threshold: float = 1.5          # alpaca_disagrees: full credit at this max(ratio_A, ratio_B)
    alpaca_mut_floor: float = 0.9          # alpaca_disagrees lowest value for a GD with >= min_support doubled mutations
    majority_coverage: float = 0.5         # minority_coverage: 0.2 at 0, 0.5 at this coverage
    minor_clone_factor: tuple[float, float] = (0.6, 0.9)  # minor_clone: [sliver of cells, >= majority_coverage]
    ambiguity_majorities: tuple[float, ...] = (0.4, 0.6)  # ambiguous_solution: majority cut-offs tried
    fixed: dict[str, float] = field(default_factory=lambda: {
        "ambiguous_solution": 0.5,
        "reversion": 0.2,
        "ploidy_path_mismatch": 0.5,
        "unplaced_events": 0.5,
        "multicluster_parent_child": 0.7,
        "fill_tie": 0.8,
        "fill_count_differs": 0.9,
        "reversion_to_diploid": 0.0,       # also reversion_to_lower_ploidy; status `failed`
    })
    reversion_min_share: float = 0.5       # a clone reverts when > this share of its cells is in samples with fewer GDs

    # majority ploidy prior (Piotr, 2026-11-12, LTX0981-Tumour1): borderline samples are moved
    # to the level most samples clearly share before the fill
    prior_enabled: bool = True
    prior_mode: str = "majority"            # "search" (2026-11-13, LTX1020) or "majority" (one-way majority move)
    prior_max_moved_fraction: float = 0.5  # search: move fewer than this fraction of all samples
    prior_max_combinations: int = 2000     # search: above this, fall back to the majority rule
    prior_band: float | None = None        # a sample is borderline within this distance of 0.5 (None: ploidy_margin)
    prior_min_clear_fraction: float = 0.5  # samples clearly at the shared level must be > this fraction of all samples
    prior_require_simpler: bool = True     # use the prior only if it gives a simpler fill than the strict targets
    prior_borderline_majority: bool = True  # also when every sample is borderline at a level and most are on one side
    # factor = lo at the edge of the band, hi just next to 0.5
    prior_raised_factor: tuple[float, float] = (0.6, 0.85)   # clonal-GD prior: raised to the doubled majority
    prior_lowered_factor: tuple[float, float] = (0.75, 0.95)  # diploid prior is stronger: lowered to the majority

    @property
    def band(self) -> float:
        return self.ploidy_margin if self.prior_band is None else self.prior_band

    def to_dict(self) -> dict:
        return asdict(self)


def load(path: Path | None) -> Params:
    """Params from a TOML file (flat keys plus a `[fixed]` table); defaults for missing keys."""
    if path is None:
        return Params()
    raw = tomllib.loads(Path(path).read_text())
    known = {f.name for f in fields(Params)}
    unknown = sorted(set(raw) - known)
    if unknown:
        raise ValueError(f"{path}: unknown parameters {unknown}")
    if "fixed" in raw:
        default_fixed = Params().fixed
        bad = sorted(set(raw["fixed"]) - set(default_fixed))
        if bad:
            raise ValueError(f"{path}: unknown [fixed] factors {bad}")
        raw["fixed"] = {**default_fixed, **raw["fixed"]}
    for k in ("ambiguity_majorities", "prior_raised_factor", "prior_lowered_factor", "minor_clone_factor"):
        if k in raw:
            raw[k] = tuple(raw[k])
    return Params(**raw)


def apply(p: Params) -> None:
    from wgd import placement, quality
    placement.MAJORITY = p.majority
    placement.MIN_COVERAGE = p.min_coverage
    placement.MINOR_CLONES = p.fill_minor_clones
    placement.MINOR_MIN_CELLS = p.minor_min_cells
    placement.ALPACA_WGD_THRESHOLD = p.alpaca_wgd_threshold
    quality.PLOIDY_MARGIN = p.ploidy_margin
    quality.MIN_SUPPORT = p.min_support
    quality.ALPACA_THRESHOLD = p.alpaca_threshold
    quality.ALPACA_MUT_FLOOR = p.alpaca_mut_floor
    quality.MAJORITY_COVERAGE = p.majority_coverage
    quality.MINOR_CLONE_FACTOR = tuple(p.minor_clone_factor)
    quality.FIXED = dict(p.fixed)
