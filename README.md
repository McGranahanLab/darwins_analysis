# darwins_analysis
Repository design to store ai drive analyses

## Analyses

- [`wgs_cn_deviation/`](wgs_cn_deviation/): TRACERx WGS single-clone tumours, deviation of fractional copy number from the nearest integer (Piotr, 2026-10-27).
- [`alpaca_post_hoc_correction/`](alpaca_post_hoc_correction/): correct ALPACA clone copy number with the agreed WGD calls, e.g. parallel GDs that ALPACA moved onto an absent parent (prototype; Piotr, 2026-11-13).
- [`WGD/`](WGD/): TRACERx whole-genome doubling. Integrated per-clone WGD table from PGDD calls, ALPACA clone CN and the CONIPHER tree. WES first, then WGS (Piotr, 2026-11-06).
- [`WGD_pipeline/`](WGD_pipeline/): bash pipeline, PGDD → WGD caller → ALPACA post-hoc correction, from tx_data, with one config per cohort: `config.wgs.toml` (WGS) and `config.wes.toml` (WES) (Piotr, 2026-11-16; WES config and rename from `WGD_WGS_pipeline/` 2026-11-24).

## Decks

- [`decks/alpaca_correction/`](decks/alpaca_correction/): how and when the post-hoc ALPACA correction is applied, with synthetic examples (Piotr, 2026-11-16).
- [`decks/wgd_pipeline/`](decks/wgd_pipeline/): rules and heuristics of the WGD pipeline, steps 00–07, with synthetic examples run through the pipeline code (Piotr, 2026-11-16).
