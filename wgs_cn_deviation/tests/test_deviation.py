"""Unit tests on hand-built frames (no tx_data needed)."""

import numpy as np
import polars as pl
import pytest

from wgscn import deviation as d


def _seg(cn_major, cn_minor, sample="S1", tumour="T1", purity=0.5, chrs=None, starts=None):
    n = len(cn_major)
    chrs = chrs or ["1"] * n
    starts = starts or [i * 10_000_000 + 1 for i in range(n)]
    return d.add_deviation(pl.DataFrame({
        "tumour_id": tumour, "patient_id": "P1", "sample": sample, "chr": chrs,
        "start": starts, "end": [s + 9_999_999 for s in starts],
        "bp": [10_000_000] * n, "cn_major": cn_major, "cn_minor": cn_minor,
        "purity": purity, "ploidy": 2.0, "coverage": 100.0,
    }))


def test_integers_have_zero_deviation():
    s = d.per_sample(_seg([1.0, 2.0, 3.0], [0.0, 1.0, 1.0]))
    assert s["dev_mean"][0] == 0
    assert s["scale_slope"][0] == pytest.approx(1) and s["scale_intercept"][0] == pytest.approx(0)


def test_half_way_is_max_deviation():
    s = d.per_sample(_seg([1.5, 2.5], [0.5, 1.5]))
    assert s["dev_mean"][0] == pytest.approx(0.5)


def test_scale_error_detected():
    true = np.array([1, 2, 3, 2, 1, 4], float)
    s = d.per_sample(_seg(list(true * 1.08), list(true * 0 + 1.0)))
    assert s["scale_slope"][0] > 1.03


def test_spearman():
    x = np.arange(10.0)
    rho, p, n = d.spearman(x, x ** 2, n_perm=2000)
    assert rho == pytest.approx(1) and p < 0.01 and n == 10
    rho, _, _ = d.spearman(x, -x, n_perm=100)
    assert rho == pytest.approx(-1)


def test_recurrent_focal_region_found():
    """A 2 Mb off-integer event shared by all tumours, each with its own segmentation."""
    rng = np.random.default_rng(0)
    frames = []
    for t in range(8):
        cuts = [c for c in set(rng.integers(1, 200, 12).tolist()) if not 100 < c < 102]
        bounds = [0] + [c * 1_000_000 for c in cuts] + [200_000_000]
        bounds = sorted(set(bounds + [100_000_000, 102_000_000]))
        starts = [b + 1 for b in bounds[:-1]]
        ends = bounds[1:]
        major = [2 + rng.normal(0, 0.05) for _ in starts]
        major = [2.45 if s == 100_000_001 else m for s, m in zip(starts, major)]
        n = len(starts)
        frames.append(d.add_deviation(pl.DataFrame({
            "tumour_id": f"T{t}", "patient_id": f"P{t}", "sample": f"S{t}", "chr": "1",
            "start": starts, "end": ends, "bp": [e - s + 1 for s, e in zip(starts, ends)],
            "cn_major": major, "cn_minor": [1.0] * n, "purity": 0.5, "ploidy": 2.0,
            "coverage": 100.0})))
    pieces = d.recurrent_pieces(d.cohort_pieces(pl.concat(frames)), n_perm=200)
    top = pieces.sort("mean_z", descending=True).row(0, named=True)
    assert (top["start"], top["end"]) == (100_000_001, 102_000_000)
    assert top["frac_tumours_z_gt_1"] == 1 and top["p_fwer"] < 0.05
    regions = d.recurrent_regions(pieces)
    assert regions.row(0, named=True)["start"] == 100_000_001


def test_region_concordance():
    rng = np.random.default_rng(1)
    base = 2 + rng.normal(0, 0.15, 50)
    a = _seg(list(base), [1.0] * 50, sample="R1")
    b = _seg(list(base + rng.normal(0, 0.01, 50)), [1.0] * 50, sample="R2")
    c = d.region_concordance(pl.concat([a, b]))
    assert c.height == 1 and c["spearman_rho"][0] > 0.9
