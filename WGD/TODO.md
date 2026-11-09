# WGD TO-DO

Deferred work. Decisions and their reasoning are in `docs/previous_analysis_logic.md`.

## Inputs (tx_data)

- [ ] Add a `conipher_mutations` source: per tumour, `<tumour>.SCoutput.CLEAN.tsv` + `treeTable.tsv` + `allTrees.txt` (WES). Needs from Piotr: the Nemo path of the WES CONIPHER output directory and a 2-line head of one SCoutput.CLEAN. (2026-11-09)
- [ ] Add a `wes_msai_segments` source: `<tumour>_copynumber_MSAI_adjusted_segments_tree_input.tsv` (WES). Needs from Piotr: the Nemo path and a 2-line head. WGS already has `wgs_consensus_segments` under a different file name, so step 00 needs a WGS path pattern. (2026-11-09)
- [ ] Then switch step 00 to read through tx_data instead of raw directories.
- [ ] Merge tx_data #6 (mock consistency, PGDD settings) and #7 (`alpaca_wgd_ratios`); retarget #7 to master after #6.

## Improvements after the first real-data run

- [x] Placement rule C: feasible-set rule (2026-11-09, LTX0025).
- [ ] Placement rule G (ploidy-only events in mutation-supported tumours): apply the same feasible-set idea.
- [ ] ALPACA ratio: review the 1.5 concerns (doc §5b): minimum `weight_sum` per allele inside `max(A,B)`, a borderline band as a quality input, and bimodality / calibration checks on real data.
- [ ] Quality: final K' components and thresholds (`N_MIN`, `M_MIN`), plus the remaining quality codes.
- [ ] Open decisions in doc §6: second GDs, triage classes, multi-cluster events with non-adjacent clusters or more than 2 clusters.

## Upstream

- [ ] PGDD `grepl` substring bug (`R/parrallel_gd_caller.R:412`): cluster 2 also matches 12 (Piotr).
- [ ] Regenerate PGDD with the new hg19 ALPACA, then hg38 (Piotr). Pass `check_mut_cpn_2_threshold = 1.25` explicitly.
