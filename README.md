# darwins_analysis
Repository design to store ai drive analyses

## Analyses

- [`wgs_cn_deviation/`](wgs_cn_deviation/): TRACERx WGS single-clone tumours, deviation of fractional copy number from the nearest integer (Piotr, 2026-10-27).
- [`alpaca_post_hoc_correction/`](alpaca_post_hoc_correction/): correct ALPACA clone copy number with the agreed WGD calls, e.g. parallel GDs that ALPACA moved onto an absent parent (prototype; Piotr, 2026-11-13).
- [`WGD/`](WGD/): TRACERx whole-genome doubling. Integrated per-clone WGD table from PGDD calls, ALPACA clone CN and the CONIPHER tree. WES first, then WGS (Piotr, 2026-11-06).

## Decks

- [`decks/alpaca_correction/`](decks/alpaca_correction/): how and when the post-hoc ALPACA correction is applied, with synthetic examples (Piotr, 2026-11-16).
