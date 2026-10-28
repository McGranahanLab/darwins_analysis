# wgs_cn_deviation

TRACERx WGS: how far is fractional copy number from the nearest integer in
**single-clone** tumours, what drives it (purity, ploidy, coverage), and is
it recurrent across or within tumours?

Requested by Piotr (2026-10-27). Built on `tx_data` canonical tables (PR #2):

| table | used for |
|---|---|
| `wgs_consensus_segments` | fractional allele-specific CN: MSAI-adjusted consensus segments; single-sample tumours use the NOT_adjusted file (`msai_adjusted` column) |
| `wgs_tree_clones` | single-clone = `n_clones_tree == 1` (clones on the ALPACA tree) |
| `wgs_battenberg_inventory` | purity, ploidy, coverage per sample |

## Measures

- Per allele, `frac = CN - round(CN)` and `dev = |frac|`: 0 = integer,
  0.5 = half-way, 0.25 = average with no integer signal. Genome summaries
  are segment-length weighted. Autosomes only by default (`--include-x`).
- **Across tumours:** per-tumour mean `dev` vs purity, ploidy and coverage
  (Spearman, permutation p). One row per tumour.
- **Purity/ploidy misfit:** `scale_slope`/`scale_intercept` fit observed
  CN on its nearest integer. A mis-scaled fit stretches every state;
  noise doesn't.
- **Within tumours:** deviation by allele × integer state, by segment
  length, and by chromosome (ratio to the sample's own mean). For
  multi-region single-clone tumours: region-to-region correlation of
  per-segment deviation.
- **Recurrent regions (segments, no bins):** the genome is cut at every
  breakpoint of every tumour, so each piece lies inside one segment per
  sample and takes that segment's deviation unchanged (z-scored within the
  sample to remove purity-driven levels). Per piece: mean z across tumours
  and the fraction of tumours with z > 1. Adjacent pieces with mean z >= 1
  are merged into regions. Genome-wide (FWER) p-values come from rotating
  each tumour's profile by an independent random offset, which keeps its
  segment structure but breaks alignment between tumours. Only pieces
  covered by >= 50% of tumours are tested (`--min-frac`).

## Run

```bash
pixi install
# build tx_data first (in the tx_data checkout, with TX_NEMO_ROOT=/nemo on Nemo):
#   pixi run python scripts/build_all.py
cd wgs_cn_deviation
pixi run python scripts/analyse.py --tx-data-dir /path/to/tx_data/data --out output/<name>
pixi run test
```

Outputs: `summary.md`, CSV tables, `fig1`-`fig4` PNGs, `run.json`
(commits + options).

## Caveats

- `minor = 0` cannot go below zero, so LOH segments have a one-sided,
  positive `frac`. Read `by_state` per allele.
- On the mock cohort (8 single-clone tumours) the analysis recovers the
  planted purity effect (rho = -1). The mock has no planted recurrent
  regions; correctly, none reaches p_fwer < 0.05 (best 0.18). A unit test
  plants a 2 Mb event shared by 8 tumours with different segmentations and
  checks it is found at its exact coordinates.
