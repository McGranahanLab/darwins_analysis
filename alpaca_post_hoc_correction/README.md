# ALPACA post-hoc correction (prototype)

Correct ALPACA clone copy number after the fact, using the agreed WGD calls from the
[WGD pipeline](../WGD/) (Piotr, 2026-11-13). Prototype on tx_data (mock) and synthetic cases.

## Why

ALPACA infers each clone's allele-specific integer copy number by parsimony over the tree and
knows nothing about genome doubling. A clone absent from every sample (an extinct ancestor)
is constrained only by its relatives. When two or more of its children doubled
**independently** (parallel GDs), the cheapest explanation for ALPACA is one shared gain on
the absent parent: the parent looks doubled (WGD ratio ~2) and the doubled children look
unchanged (~1). Seen in LTX1348 (clone1 1.78; clone5 + clone6), LTX0552-T2 (clone1 2.14; six
branches), LTX0501 (clone1 1.95) and LTX0929 (clone2 2.07). The WGD pipeline currently works
around it when scoring (`alpaca_parent` flag); this project fixes the copy number itself.

## Rule

A clone is corrected when it is **absent from every sample**, carries **no GD** in the calls,
and a **GD-carrying clone is reached from it through absent, GD-free clones only**: a child with a
GD, or a chain of missing clones above one (Piotr, 2026-11-13: if clone1 is absent and its child
clone2 is absent too, and clone2's children doubled in parallel, both clone2 and clone1 are
corrected, and nothing else; an observed clone stops the chain). Clones are processed deepest first, so corrected values
feed the next level up. For every segment and allele, the new copy number `x` minimises

    sum over children c  |CN_c - 2^g_c * x|   +   |x - CN_parent|

with `g_c` the GDs on child `c` and `CN_parent` the parent's copy number (1 for `diploid`).
Ties keep ALPACA's value if it is among the minimisers, else the smallest. Clones with cells in
any sample are observed and never changed. WGD ratios are then recomputed exactly as ALPACA's
`analysis.calculate_wgd_ratios` (checked against tx_data's `alpaca_wgd_ratios`).

Example (test `test_parallel_gds_absent_parent_halved`, LTX1348-like): ALPACA clone1 2|2, 2|2,
2|2, 4|2 → corrected 1|1, 1|1, 1|1, 2|1. clone1's ratio_A goes 2.67 → 1.33 and clone5's 1.08 → 2.17,
so the ratios now agree with the GD calls.

## Run

```bash
pixi install
pixi run correct --tx-data-dir ../../tx_data/data --out output/mock \
    [--wgd-calls ../WGD/output/<name>/05_table/wgd_calls.tsv] [--min-quality 0.3]
pixi run test       # needs tx_data next to darwins_analysis, or TX_DATA_REPO
```

Inputs: tx_data `alpaca`, `clone_proportions`, `alpaca_wgd_ratios` (tree). WGD calls:
`--wgd-calls` (WGD step 05 `wgd_calls.tsv`), else tx_data `wgd_calls` (still the placeholder).
`--min-quality`: only tumours whose call quality reaches it are corrected; the rest are listed
in `skipped.tsv`.

Outputs (`<out>/`): `alpaca_corrected.parquet` (tx_data `alpaca` schema, derived columns
recomputed, plus `pred_CN_A_original`, `pred_CN_B_original`, `corrected`),
`wgd_ratios_corrected.tsv` (ALPACA's 12 columns), `ratio_comparison.tsv`, `corrections.tsv`
(one row per corrected clone: children with GDs, bp fraction changed, ratios before/after),
`skipped.tsv`, `flags.tsv`, `run.json`.

## Flags (not fixed)

Disagreements left after the correction are flagged only (Piotr, 2026-11-13: they should never
happen): `gd_without_alpaca_doubling` (a clone with a called GD whose ALPACA ratio vs its parent,
max(A, B), is below 1.5) and `alpaca_reversion` (a clone below a GD-carrying clone whose
length-weighted median total CN is ≤ 2). On the tx_data mock, with the placeholder `wgd_calls`,
two clones are flagged; that only shows the mechanics.

## Status

- tx_data mock (9 tumours): no clone meets the rule, so nothing changes. The ratio port matches
  ALPACA to 1e-9 on all 34 clones.
- Real cases need the real ALPACA output (tx_data `alpaca` is a mock) and the released WGD calls.

## Decisions (Piotr, 2026-11-13)

- Only missing clones are corrected (single clones or chains of them); observed clones never, at
  least until live cases are checked.
- ALPACA vs call disagreements on observed clones and reversions: flag only.
- The corrected tables go to tx_data as new sources (`alpaca_corrected`,
  `alpaca_wgd_ratios_corrected`).
