# WGD TO-DO

Deferred work. Decisions and their reasoning are in `docs/previous_analysis_logic.md`.

## Inputs (tx_data)

- [ ] Add a `conipher_mutations` source: per tumour, `<tumour>.SCoutput.CLEAN.tsv` + `treeTable.tsv` + `allTrees.txt` (WES). Needs from Piotr: the Nemo path of the WES CONIPHER output directory and a 2-line head of one SCoutput.CLEAN. (2026-11-09)
- [ ] Add a `wes_msai_segments` source: `<tumour>_copynumber_MSAI_adjusted_segments_tree_input.tsv` (WES). Needs from Piotr: the Nemo path and a 2-line head. WGS already has `wgs_consensus_segments` under a different file name, so step 00 needs a WGS path pattern. (2026-11-09)
- [ ] Then switch step 00 to read through tx_data instead of raw directories.
- [ ] Merge tx_data #6 (mock consistency, PGDD settings) and #7 (`alpaca_wgd_ratios`); retarget #7 to master after #6.

## Improvements after the first real-data run

- [x] Placement rules C and G: replaced by the tumour-level ploidy fill (2026-11-09, LTX0025, LTX0073).
- [ ] Ploidy fill: review the 50% sample-level cut-off and the tie-breaks (spill, then ALPACA score) on the real cohort; check `fill_tie` events.
- [ ] ALPACA ratio: review the 1.5 concerns (doc §5b): minimum `weight_sum` per allele inside `max(A,B)`, a borderline band as a quality input, and bimodality / calibration checks on real data.
- [ ] Highlight cases for review (Piotr, 2026-11-10). The pipeline should list, per tumour, with the reason:
  - a) **failed QC**: excluded by manual QC (`FINALIZED`) or truncal CCF ≤ 90. Today these are silently dropped in step 04, so report them instead (tumour, which check failed);
  - b) **low quality**: `quality = Low`, with its `quality_code`;
  - c) **unlikely solutions**, e.g. reversion to diploid (step 04 flag `reversion`). Also candidates: `fill` incomplete (`missing:`/`over:`), `fill_tie`, `no_alpaca_tree`, `fixed_event_not_on_tree`, the path-vs-ploidy mismatch (H), a fill GD count differing from PGDD's ploidy-only event count, and multi-cluster parent–child events.
  Suggested output: `05_table/review.tsv` (tumour_id, category, reason) plus counts in a run summary.
- [ ] Quality: final K' components and thresholds (`N_MIN`, `M_MIN`), plus the remaining quality codes.
- [ ] Open decisions in doc §6: second GDs, triage classes, multi-cluster events with non-adjacent clusters or more than 2 clusters.

## Upstream

- [ ] PGDD `grepl` substring bug (`R/parrallel_gd_caller.R:412`): cluster 2 also matches 12 (Piotr).
- [ ] Regenerate PGDD with the new hg19 ALPACA, then hg38 (Piotr). Pass `check_mut_cpn_2_threshold = 1.25` explicitly.
