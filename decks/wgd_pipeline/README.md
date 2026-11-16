# Deck: WGD pipeline rules and heuristics

The rules and heuristics of [`WGD/`](../../WGD/), steps 00–07, for the lab (Piotr, 2026-11-16).
17 slides, speaker notes on each: `wgd_pipeline.pptx`.

| # | Slide |
|---|---|
| 1–2 | Title; the eight steps and their inputs |
| 3 | Step 00: GDs per sample from rounded major CN (computed by `compute_num_gds`) |
| 4 | Steps 01–02: PGDD event kinds, bootstrap-stable calls, pruning to the ALPACA tree |
| 5 | QC: upstream QC and CCF QC, `qc_failed` / `no_alpaca_tree` rows (PR #6) |
| 6 | Triage classes, class and sub_class |
| 7 | Fixed GDs: clonal → MRCA, mutation-supported → its clone, parent–child pair → child |
| 8 | Why not per event: the old rule over-doubles a sample, the fill doesn't |
| 9–10 | The ploidy fill and its tie-breaks |
| 11 | Relaxed and minor-clone passes |
| 12 | The ploidy prior for borderline samples |
| 13 | Failures and ALPACA checks (reversion detector, H, unplaced events) |
| 14–15 | Quality factors (values read from `wgd.params`); per-event quality and `--reject-below` |
| 16 | Outputs and the review list |
| 17 | Tuning, safeguards, open questions |

## Where the numbers come from

`examples.py` builds synthetic tumours (`SYN0001-*`) and runs them through the pipeline's own step 04
and 05 code (`wgd.placement.place_events`, `wgd.final_table.build_table`); slide 3 runs
`wgd.pgdd_input.compute_num_gds`; slide 14 reads the `wgd.params` defaults. `tests/` checks the outcomes
on the slides against the code and that no real tumour id appears. Shared slide primitives:
[`../deckkit.py`](../deckkit.py).

## Rebuild

```bash
pixi install
pixi run deck      # wgd_pipeline.pptx
pixi run test
pixi run preview   # rough PNGs + text-overflow report in preview/ (no LibreOffice needed)
```

Needs `../../WGD` (imported, not installed). The preview is approximate; the deck has not been opened
in PowerPoint yet.
