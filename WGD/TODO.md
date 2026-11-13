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
- [x] Quality: continuous per-tumour score (doc K'', `wgd/quality.py`, 2026-11-10).
- [ ] Quality: calibrate weights on the cohort (defaults accepted for now, 2026-11-10). Per-event `borderline_ploidy` is now level-specific (2026-11-12, LTX0693); the tumour-level factor still uses every level. On the real cases, `alpaca_disagrees` fires for most subclonal GDs (ratios 1.0–1.3), so check whether ALPACA ratios are informative for subclonal GDs before keeping that factor at full weight.
- [ ] Each new edge case: copy its files into tests/cases/<tumour>/ (cp_table, ratio_scores, pgdd rows), add the agreed answer to agreed.toml.
- [ ] Majority ploidy prior (2026-11-12, LTX0981): count `ploidy_prior` applied / not_simpler on the real cohort; check that lowering 2→1 should use the diploid (stronger) factor; calibrate `prior_*_factor` (defaults accepted for now).
- [ ] Quality: score the *amount* of doubling, not only the genome fraction past the cut-off (Piotr, 2026-11-12). Today a sample counts as doubled when > 50% of the genome has major CN rounding to >= 2. Rounding is half-to-even (`pgdd_input.compute_num_gds`), so the real cut-off is raw major CN >= 1.5. A sample at major CN 1.6 genome-wide therefore scores the same as one at 2.0, though the second is more believable. Option: a graded factor on the length-weighted mean major CN (or its distance from the level's target) per sample. Step 00 already computes `conf_gd1/2/3` (`_level_confidence`: length-weighted closeness of (A|B) to (2|2), (4|4), (8|8)), but nothing downstream uses it. Check whether it fits (it uses both alleles, so LOH lowers it) or whether a major-allele-only version is better. Calibrate on the cohort.
- [ ] Open decisions in doc §6: second GDs, triage classes, multi-cluster events with non-adjacent clusters or more than 2 clusters.

## Post-hoc ALPACA correction

- [ ] When our GD calls contradict ALPACA's parsimony, correct ALPACA's clone CN afterwards (Piotr, 2026-11-12). ALPACA doesn't know PGDD: if an absent clone has several doubled children, it moves the shared gain onto that parent (one event instead of parallel GDs). Once we know the parent has no GD, its CN should be roughly halved, and the children's ratios then agree. Cases: LTX1348 (clone1 1.78; clone5 + clone6), LTX0552-T2 (clone1 2.14; six branches), LTX0501 (clone1 1.95; clone4 + clone5 + clone6). Until then, quality uses the parent's ratio for such GDs (`alpaca_parent` flag). This is a separate step after placement, not part of the scoring.

## Micropipeline

- [ ] Wire the modules into one micropipeline (Piotr, 2026-11-13): WGD steps 00–07 → alpaca_post_hoc_correction (darwins_analysis#3) → tx_data build (`alpaca_corrected`, tx_data#8). Every input and output location is a parameter; nothing hardcoded.

## Upstream: PGDD

- [ ] PGDD `grepl` substring bug (`R/parrallel_gd_caller.R:412`): cluster 2 also matches 12 (Piotr).
- [ ] Mutation-supported events include samples that contradict them (Piotr, 2026-11-12: a PGDD issue, not ours; no change to our method for now). **LTX0969-T1:** GD1 is on cluster 2 in BR_T1.R1 + SU_T1.R1, but in SU_T1.R1 only 8/106 cluster-2 mutations are doubled (7.5%; discovery fails, the looser check step passes), against 23/36 (64%) in BR_T1.R1. SU_T1.R1 has 1 GD (`nMaj3` 0.20), so its GD precedes cluster 2, i.e. it's clonal. With GD1 fixed on clone2, a clonal GD would put SU_T1.R1 at 2, so the fill gives clone2 + clone3 + clone4 ×2 (4 GDs). Moving GD1 below clone2 (BR_T1.R1 only) gives clone1 + clone3 (or clone8) + clone4: 3 GDs, every sample at its PGDD level. This also fits ALPACA (clone1 1.97, clone4 1.63, clone3 1.29, clone2 1.01), truncal mutations 71–79% doubled in all samples, and cluster 4 7/10 doubled in SU_T1.R2. Cohort: 13 events in 13 tumours have a sample with ≥ 10 cluster mutations present but < 25% doubled (LTX0063, 0552-T2, 0572 (6 samples), 0731, 0857, 0944, 0969, 0979, 0982, 1000, 1253, 1271, 1281). None of the regression cases are among them. Options: fix in PGDD (require the discovery criterion per sample, not just the check), or, on our side, drop contradicting samples from the event and move the GD down the cluster's subtree (proposed 2026-11-12, not implemented).
- [ ] Regenerate PGDD with the new hg19 ALPACA, then hg38 (Piotr). Pass `check_mut_cpn_2_threshold = 1.25` explicitly.
