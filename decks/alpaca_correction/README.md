# Deck: post-hoc ALPACA correction

How and when [`alpaca_post_hoc_correction`](../../alpaca_post_hoc_correction/) re-infers ALPACA's
copy number for absent ancestors after parallel genome doublings, with synthetic examples
(Piotr, 2026-11-16). 12 slides, speaker notes on each: `alpaca_correction.pptx`.

| # | Slide |
|---|---|
| 1 | Title |
| 2 | The problem: ALPACA can't see a GD, so parallel GDs become one gain on the absent parent |
| 3 | When: tumour gates (calls, tree, `--min-quality`) and the three clone conditions |
| 4 | Six synthetic trees: which clones qualify (A, C, C2 corrected; D1–D3 not) |
| 5 | How: the per-segment, per-allele cost and the tie rule |
| 6–7 | Example A (parallel GDs): cost table, corrected CN, ratios before/after |
| 8 | Example B: an undoubled child keeps ALPACA's value |
| 9 | Example C: chains corrected deepest first, and why the order matters |
| 10 | What is never changed; what is only flagged |
| 11 | Outputs and tx_data sources |
| 12 | Status and next steps |

## Where the numbers come from

`examples.py` builds the synthetic tumours (`SYN-*`, 4 segments) and runs them through the real
`apc.correct.run`; the cost tables use the same objective as `correct.best_cn`. `make_deck.py` draws
everything from those results (native, editable shapes, tables and a chart). `tests/` checks the
outcomes and numbers on the slides against the code and that no real tumour id appears.

## Rebuild

```bash
pixi install
pixi run deck      # alpaca_correction.pptx
pixi run test
pixi run preview   # rough PNGs + text-overflow report in preview/ (no LibreOffice needed)
```

Needs `../../alpaca_post_hoc_correction` (imported, not installed). The preview is approximate
(DejaVu fonts scaled to Calibri/Cambria widths); the deck has not been opened in PowerPoint yet.
