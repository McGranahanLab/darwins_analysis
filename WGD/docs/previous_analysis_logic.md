# Logic of the previous WGD consolidation (`templates/previous_analysis.ipynb`)

How the previous run turned PGDD calls, ALPACA output and QC tables into
per-clone WGD calls, written out as a pipeline. Cell numbers `[n]` refer to
the notebook (58 cells, no saved outputs). Where the code and its comments
disagree, this document follows the code and flags the difference.

Written by Darwin, 2026-11-06, for the WGD rebuild; revised 2026-11-07 after
Piotr's review. Section 5 lists what probably should **not** be carried over.
Section 6 records the decisions. Section 7 lists open checks and reminders.

---

## 1. Overview

```
PGDD run ─┐
          ├─ A. tumour QC (manual QC table + truncal CCF)
QC tables ┘
          └─ B. triage tumours (GDs_per_tumour)
               ├─ no_gd ..................... resolved: GD = 0 everywhere
               ├─ clonal_gd_only ............ resolved: num_clonal_gds on MRCA
               ├─ ploidy-only tumours ....... C. ALPACA ratio_A assignment
               └─ mutation-supported tumours  D. manual review queue
          E. multi-cluster events -> one clone (parent–child pairs only)
          F. clonal events -> MRCA
          G. ploidy-only events -> clone, from tree position (single-sample events only)
          H. QC: GDs per root→leaf path vs max sample num_gds
          I. QC: reversion (descendant clone with median CN_A = 1)
          J. draft per-clone table (unused) ── manual reconciliation in wgd_viewer ──
          K. quality from supporting mutation counts -> final table (from viewer output)
```

The final table (K) has the same shape as the placeholder `wgd_calls` in
tx_data: `tumour_id, clone, parent, number_of_gds, total_gds, quality, status, class`.

## 2. Inputs

| Input | Path in notebook | Used for |
|---|---|---|
| PGDD outputs | `CN-CCF/pgdd/wes/{release}/pgdd_output/{run_name}/` with `run_name = hypergeometric_relaxed`: `pgdd_GD_events`, `_GDs_per_region`, `_GDs_per_tumour`, `_mut_counts`, `pgdd_cohort_gd_summary` | Everything |
| PGDD input | `.../pgdd_input/{run_name}/pgdd_input_all_tumours.tsv` | Loaded, not used |
| `pgdd_cohort_gd_summary.tsv` | as above | **Used in H.** It supplies the per-sample `num_gds` (key `sample_name_hash`). The same values exist as `num_gds` in `GDs_per_region`. |
| Pipeline PGDD results (2025-08-09) | `tracerx/_RELEASE/release_tx842/_aggregate/pgddetect/20250809_cohort_*` | Loaded and QC-filtered for comparison only. No logic depends on them. |
| ALPACA per-tumour WGD scores | `tracerx800/output/tx842/snp_filt_0_cn_rec_FFT/tumour_results/<tumour>/wgd_ratio_scores.csv`: `tumour_id, clone, parent, ratio_A` | **The tree** (clone→parent) and a per-clone WGD score (`ratio_A`). Written by `alpaca wgd` (McGranahanLab/ALPACA-model, `alpaca/analysis.py:120 calculate_wgd_ratios`): `ratio_A` is the segment-length-weighted mean of the clone's allele-A CN divided by its parent's, over segments where the parent's CN > 0. The MRCA's parent is `diploid` (CN 1+1). The current version also writes `ratio_B`, `ratio_total` and segment counts/weights. Being added to tx_data. |
| ALPACA combined CN | `.../cohort_results/all_tumours_combined.csv` (`pred_CN_A`, `pred_CN_B`) | I (reversion) and the J template (clone list) |
| ALPACA clone proportions | `tracerx800/input/tx842/snp_filt_0_cn_rec_FFT/<tumour>/cp_table.csv` (clone × sample) | Which clones are present in an event's samples (C, G) |
| Manual QC table | `_aggregate/20250729_tx842_manual_qc_complete_post_refphase.tsv` | A |
| CONIPHER `clusterInfo.txt` | `_RELEASE/release_tx842/<patient>/mutation_trees/tumour_<n>/clusterInfo.txt` | A (truncal CCF) |
| Problematic events | `CN-CCF/pgdd/_assets/prototype_events.tsv` (`tumour_max_problematic_level`) | D (manual cross-check only) |
| wgd_viewer output | `CN-CCF/pgdd/wes/release_tx842/latest_viewer_output/clone_gd_calls_latest.tsv` | K (result of manual reconciliation) |
| Queue files (written) | `CN-CCF/pgdd/queue_files/pgdd_wes_{no_gd,clonal_gd_only,alpaca_called,mut_supported}_queue.tsv` | Status per tumour. The queue name becomes `class` in the final table. |

The PGDD outputs differ from the current run in two ways:
- **`gd_p_value`** (a hypergeometric test, `run_name = hypergeometric_relaxed`) is a column of `mut_counts` in the old run. **We no longer use the hypergeometric test.** The current PGDD (`cb0437b`) uses thresholds plus bootstrap stability, and writes only stable calls to the `pgdd_*` files; unstable calls go to `unstable_calls_pgdd_*`. Step 01 reads only the `pgdd_*` files, so every call we see is already stable.
- **`clusters`** held strings like `clone3` (cell `[6]` splits on `"clone"`). The current output writes `3`.

## 3. Stage-by-stage logic

### A. Tumour QC `[3] get_qc_df, additional_ccf_qc; [5]`
1. **Manual QC.** From the manual QC table, keep rows with `selected_solution == True`. Build the tumour id as `patient + "-Tumour" + tumour_id`. A tumour passes only if **all** its rows have `FINALIZED == True`.
2. **Truncal CCF QC.** For each passing tumour, read CONIPHER `clusterInfo.txt` and keep clusters with `treeClust & truncal`. Take the **lowest** `meanCCF` among truncal clusters. The tumour passes if that value is above **0.9**. `meanCCF` is in percent (0 to slightly above 100), so the port uses **90** (Piotr, 2026-11-09: the intended cut-off is 90%; the fraction conversion happened elsewhere in the old run). Missing `meanCCF` only occurs on clusters that are off the tree.
3. All PGDD tables are subset to tumours that pass both checks. The cohort summary is subset to their samples.

### B. Tumour triage `[11], [14]`, from `GDs_per_tumour`
| Class | Rule | Outcome |
|---|---|---|
| `no_gd` | `First_GD == "No GD"` | resolved. Every clone has GD = 0. |
| `clonal_gd_only` | `First_GD == "Clonal"` and `num_subclonal_gds == 0` | resolved. `num_clonal_gds` goes on the MRCA. |
| ploidy-only (`alpaca_called`) | the rest, where **no** event in the tumour has `is_subclonal_mutation_supported` | C |
| `mut_supported` | the rest, where **at least one** event has mutation support | D, manual |

Note that a tumour whose `First_GD` is subclonal but which also has a clonal second GD falls into the "rest" group.

### C. Ploidy-only tumours: assignment by ALPACA ratio `[14], [15]`
For each subclonal event, in file order:
1. Find the clones present in the event's samples: clones with `cp_table` proportion > 0 in any of those samples (summed across them). Add **all ancestors** of those clones, since extinct ancestors can carry the event.
2. Candidates are the ALPACA `wgd_ratio_scores` rows for those clones, **excluding the MRCA** (`parent != "diploid"`) and excluding clones already given to an earlier event in the same tumour.
3. Pick the candidate with the **highest `ratio_A`**. (New: we use `max(ratio_A, ratio_B)`, see decision 2.) Write `clusters = "alpaca_<clone>"` and store `alpaca_ratio_A`.

The tumour is `resolved` if **every** event has `ratio_A > 1.5`. Otherwise it is `unresolved` and gets checked in wgd_viewer.

### D. Mutation-supported tumours `[8], [16], [18]–[22]`
- Every such tumour goes into the queue as `unresolved` and is reviewed manually.
- **Obsolete in the old notebook:** the review memo showed, per event, the 2×2 table behind the hypergeometric `gd_p_value` and called an event `significant` if `gd_p_value < 0.05`. We don't run that test any more.
- **Now:** PGDD's bootstrap decides which calls are kept. Per (tumour × sample × cluster) it resamples mutations (`n_boot = 500`) and keeps a call if its stability passes (`stability_threshold = 0.7`, binomial test, BH-adjusted, α = 0.05; see `pgdd_bootstrap_discovery_settings.tsv`). Unstable calls are written to `unstable_calls_pgdd_*` and never reach us. So there is no significance filter to apply in our pipeline; `stability` / `q_value` from `pgdd_bootstrap_discovery_summary.tsv` are available as evidence (e.g. for quality, K).
- Tumours are cross-checked against `prototype_events.tsv` (previously flagged events) and a hard-coded list of 12 already-reviewed tumours. This is bookkeeping only.
- `[21]` describes an idea that was never implemented: apply step C to the ploidy-only events inside mutation-supported tumours. There, a low `ratio_A` would mean the *timing* is uncertain, not the event's existence. `[22]` gives target examples: LTX0019-T1 should go to clone3, and LTX0064-T1 to clones 32 and 1.

### E. Multi-cluster events `[24]`
- Events with exactly 2 clusters: if the two clusters are a direct **parent–child** pair on the ALPACA tree, the event is assigned to one of them (`clusters = "multicluster_<clone>"`).
  - ⚠️ The comment says "move the event to the parent", but the code selects the row where `clone` is one cluster and `parent` is the other. That row's `clone` is the **child**, so the code assigns the event to the child.
  - **Decision (Piotr, 2026-11-09): the child.** If the event were placed on the parent, the doubled mutations of the child could only be explained by a second GD in the child. A parent–child pair therefore means something went wrong: either (a) the two clusters should be one (wrong split), or (b) the child really has a second GD, which should then show in its ploidy. Such events get a quality penalty (flag `multicluster_parent_child`, K').
- Events with more than 2 clusters are skipped with a warning.
- Pairs that aren't parent–child are left unchanged. The comment says they should be "flagged", but no flag is set.

### F. Clonal events `[25]`
Every event with `is_clonal` goes on the MRCA, the ALPACA clone with `parent == "diploid"`, as `clusters = "clonal_<mrca>"`.

### G. Ploidy-only events: placement by tree position `[26]`
Applies to subclonal events with no clusters that were not already resolved by C. **In practice this means ploidy-only events inside mutation-supported tumours only:** C has already written `alpaca_<clone>` (or `alpaca_None`) into every subclonal event of the ploidy-only tumours, so those events no longer have empty `clusters` and never reach G.

**Idea in one sentence:** put the event on the most ancestral clone that could explain the doubled sample without colliding with the tumour's other (already placed) GDs. A clone collides if it is another event's clone, its ancestor (the other event would then be redundant) or its descendant (this would be a second GD, which is out of scope here).
1. Skip the tumour if it has any second GD (`num_second_gd > 0`).
2. Skip the event if any other event in the tumour is ploidy-based (unplaced) or multi-cluster.
3. Collect the clones the event can't go on. These are the clones of the other events in the tumour, plus their ancestors and descendants.
4. Only events seen in **one sample** are handled. For these, ALPACA ratios are deliberately not used: the MRCA is usually absent from that sample, and ALPACA then resolves the region parsimoniously as a clonal GD. Multi-sample events are skipped.
5. Take the clones present in the sample plus their ancestors (not `diploid`). Split them into lineages, each ordered from leaf to root.
6. Walk each lineage from the leaf toward the root. A clone is a candidate if it is allowed and none of its descendants are disallowed. The leaf is a candidate as long as it is allowed.
7. In each lineage, take the **last** candidate, which is the most ancestral allowed clone.
   - With one lineage, use that clone.
   - With several lineages, use their shared top clone if they all end on the same one. Otherwise leave the event unresolved.
   - With no candidates at all, mark the event `false_call` (a suspected PGDD error).
   - Write `clusters = "ploidy_<clone>"`.

What the comment intended, but the code doesn't do, is a third condition: the chosen clone's whole subtree must sit in the doubled sample (section 5, item 5). Without it the walk only stops at collisions, so it tends to climb as high as possible. In a mutation-supported tumour the other event's ancestors include the MRCA, so the MRCA is always excluded and the climb stops below it.

### H. QC: tree-path GD count vs ploidy `[27]`
For each tumour that has subclonal events, all of them placed:
- For every leaf, count the GD-carrying clones on the leaf→root path.
- Compare the **maximum** count over leaves with the **maximum** per-sample `num_gds` (from the cohort summary).
- If they differ, flag the tumour for manual review.

This came out of LTX1167-T1. There, a mutation-supported GD on clone2 already explained the ploidy, so PGDD's extra ploidy-only event was spurious. Tumours with unplaced or `false_call` events are skipped.

### I. QC: reversion `[31], [35]`
For the subclonal events of mutation-supported tumours:
- Take the event's clone and all its descendants.
- If any descendant has a **median `pred_CN_A` of 1** across its segments, tag the event `reversion_…`. This suggests the clone has returned to a diploid-like state after the WGD.

`[35]` starts a second check for ploidy-only calls: compare the lowest `ratio_A` among the MRCA's children with the MRCA's own `ratio_A`. It wasn't finished, and it depends on a leftover `mrca` variable.

### J. Draft per-clone output table, then manual review in wgd_viewer `[36]–[43]`
"Template" here only means the **empty skeleton of the final table**: one row per (tumour, ALPACA clone), with every clone set to `GD = 0` and `Quality = NaN`, which the notebook then fills in where it can.
- Add the QC flags: `tx_QC` (step A1), `ccf_QC` (step A2) and `QC_pass = tx_QC & ccf_QC`.
- Fill in the two automatic classes:
  - `clonal_gd_only`: put `num_clonal_gds` on the MRCA, `Quality = 1`.
  - `no_gd`: all clones stay 0, `Quality = 1`.
- **This table is a dead end in the notebook:** its export (`june_reconciliation.tsv`, `[43]`) is commented out, and nothing later reads it.

The real hand-over is the manual step: every tumour with subclonal events is reconciled by hand in **wgd_viewer** (the lab's interactive app), working from the queue files and the events table. The viewer's "generate output" button writes `clone_gd_calls_latest.tsv`, which is copied to Nemo by hand (`[48]`). Step K starts from that file, not from the draft table. How the viewer builds its own per-clone rows isn't visible in the notebook.

### K. Quality and final table `[44]–[57]`
1. **Supporting mutation count.** For each mutation-supported event, look at the `mut_counts` rows for its cluster that pass `discover_pass_all` or `check_pass_all`. Take the **maximum `num_cn2`** across those rows. For multi-cluster events, the code uses the single resolved cluster from E.
2. Merge the counts onto the viewer output by `(tumour_id, clone_id)`.
3. If any event in a tumour has a supporting count at or below the **5th percentile** across events, set that tumour's quality to 0.5. The viewer's own quality is kept where it was set, and missing values become 1.
4. Tumour quality is the **minimum** over its clones, then mapped to `Low` (< 1) or `High`.
5. Output columns: `tumour_id, clone_id, number_of_gds, quality, status, class (= viewer queue), total_gds (sum per tumour), parent` (parent taken from ALPACA `wgd_ratio_scores`).

`[48]` mentions also computing the "Bielski fraction" for quality. This wasn't implemented.

#### K''. Implemented quality score (2026-11-10, Piotr: continuous, per tumour, no GD order)

Step 04 computes `quality` in [0, 1] per tumour, as the product of factors (`wgd/quality.py`). Every factor below 1 is listed in `quality_code` with the value that triggered it. Per-event flags stay in `events_placed.tsv`.
- **Graded, 0.5 → 1:**
  - `borderline_ploidy`: the smallest distance of any sample's deciding `genome_frac_nMaj*` from 0.5; full credit at ≥ 0.10.
  - `low_mut_support`: the smallest supporting count of a mutation-supported GD; full at ≥ 10.
  - `alpaca_disagrees`: the smallest `max(ratio_A, ratio_B)` among GD clones; 0.5 at 1, full at ≥ 1.5.
- **Fixed:**
  - `ambiguous_solution` 0.5: the fill answer changes when a borderline sample is moved across 0.5 (each alone, then all together) or the clone-majority cut-off is 0.4/0.6. The alternatives are listed.
  - `reversion` 0.2: a clone below a doubled clone has median ALPACA total CN ≤ 2.
  - `ploidy_path_mismatch` 0.5, `unplaced_events` 0.5, `multicluster_parent_child` 0.7, `fill_tie` 0.8, `fill_count_differs` 0.9.
- **The old K rule** (cohort 5th percentile of supporting counts) is dropped.
- **Per-event score** (2026-11-11): `worst_event_quality` / `worst_event_quality_code` on each GD-carrying clone in step 05 (the worst event when a clone carries several), and per event in `events_placed.tsv`. These are the factors that belong to one GD: ploidy margin over its samples, its mutation support, its clone's ALPACA ratio, minority coverage, tie, multi-cluster, whether its clone's GD count changes in an alternative solution, and reversion charged to the nearest doubled ancestor. Path mismatch, unplaced events and fill-count difference are tumour-level only. A clone with several GDs shows its worst event.
- **Per-event table and rejection** (2026-11-11): step 05 writes `wgd_events.tsv`, one row per GD with id `tumour:clone:k`. The tumour score is computed as `Q.aggregate`: for each factor code, the worst value over the tumour's events and its tumour-level factors (borderline ploidy over all samples, path mismatch, unplaced events, fill-count difference, and ambiguity when no event carries it), then the product. This equals the earlier definition, but can be recomputed after `--reject-below Q` drops GDs with event quality < Q. One cut-off is used for the whole cohort. Mutation-supported GDs are never rejected unless `--reject-mutation-supported` is given (Piotr, 2026-11-11: PGDD's bootstrap already triaged them). Each rejection is listed in `rejected_events_below_<Q>.tsv`. Rejected GDs are listed (`rejected = true`, `n_rejected_events`) and not penalised: the remaining calls are scored as they stand. The samples they explained are then under-explained, which a reviewer should know about.
- **Status:** defaults accepted as a starting point (Piotr, 2026-11-10). Weights and thresholds are tuned later on the cohort. Open questions: whether `alpaca_disagrees` is informative for subclonal GDs (real subclonal GD clones show ratios of 1.0–1.3), and whether `borderline_ploidy` should look only at the first GD level (LTX0025 and LTX0073 score low because of one sample about 0.03–0.04 from 0.5 on a deeper level).

#### K'. Proposed replacement quality (draft, 2026-11-07; superseded by K'')
Piotr: the metric needs improving, with three components. Proposal: compute each component **per event**, keep them as columns, and derive the label from fixed thresholds (not cohort percentiles).

| Component | Per-event value | Flag when | Source |
|---|---|---|---|
| **Mutation support** | `n_support` = `num_cn2` of the supporting cluster (max over the event's regions); also `stability` / `q_value` from the bootstrap summary | `n_support < N_MIN` (e.g. 10, TBD). Ploidy-only events have no mutation evidence and are flagged `no_mut_support` | `pgdd_mut_counts`, `pgdd_bootstrap_discovery_summary` |
| **Borderline ploidy** | `ploidy_margin` = min over the event's regions of (fraction of genome with major CN ≥ 2) − 0.5. For a second GD the fraction uses major CN ≥ 3. This is the "Bielski fraction" of `[48]` | `ploidy_margin < M_MIN` (e.g. 0.05–0.10, so 55% doubled is borderline). Arguably also check the regions *without* the event: 45% is a borderline "no". | Not in the PGDD outputs: PGDD takes `num_gds` as an input. Source (2026-11-09): step 00 computes it from the MSAI-adjusted segments (`frac_maj_ge2/3/4`), and step 03 writes it per sample as `genome_frac_nMaj2/3/4` in `pgdd_cohort_gd_summary.tsv`. |
| **Parent–child multi-cluster** | `multicluster_parent_child` = the event's PGDD clusters include a direct parent–child pair (placed on the child, E) | always (wrong cluster split, or an unseen second GD in the child) | `pgdd_GD_events.clusters` + ALPACA tree |
| **GD order** | `gd_order` = 1 for a first GD, 2 for a second, 3 for a third, counted along the root→leaf lineage | `gd_order ≥ 2`, stronger penalty for 3 | Our placement (04) |

Quality codes: every quality below High carries `quality_code`, a `;`-joined list of reasons. Current codes:
- `single_clone_num_gds_discordant`: a single-clone tumour whose samples disagree on `num_gds`. The mode is used and quality is **Low** (Piotr, 2026-11-09).
- `num_gds_mode_tie`: no unique mode; the lowest tied value is used.
- Planned: `multicluster_parent_child`, `low_mut_support`, `no_mut_support`, `borderline_ploidy`, `gd_order_2` / `gd_order_3`, `alpaca_ratio_borderline`.

Combining: event quality = `High` with no flags, `Medium` with one, `Low` with two or more (or with `no_mut_support` plus `borderline_ploidy`). Tumour quality = the worst event, but the per-event flags are kept so users can filter on the reason. Thresholds `N_MIN` and `M_MIN` are open (section 7, reminder 4).

## 4. How it maps to the current WGD steps

| Previous | Current |
|---|---|
| A. QC | not yet. It needs the QC table and CONIPHER `clusterInfo`, or tx_data equivalents. |
| (implicit) | **02** pruning: off-tree clusters and tumours without a tree. The notebook never handles off-tree clusters explicitly. `get_descendants` on an off-tree clone just returns `[]`. |
| B. triage | 04, derived from `events_pruned.support` rather than PGDD's per-tumour flags |
| E, F, G, C | 04 placement |
| H, I | 05 checks against ALPACA |
| J, K | 06 integrated table |

### Port notes (steps 04–06, 2026-11-09)

The notebook's automatic rules are ported as is in `wgd/placement.py` (step 04) so the new table can be checked against the previous release (step 06). Found while porting:
- **Rows** come from wgd_viewer in the old flow. Here they are generated: one row per ALPACA clone (`wgd_ratio_scores`) of every QC-passed PGDD tumour. `status` is the notebook's automatic queue status: `resolved` for no_gd, clonal_gd_only, and alpaca_called with every ratio > 1.5; `unresolved` otherwise. The manual statuses and the viewer's own quality can't be reproduced.
- **E with old cluster labels.** The notebook builds `f"clone{c}"` from `clusters`. If the old run wrote `clone3,clone6` (cell [6] suggests it did), that gives `cloneclone3` and the multi-cluster rule never fires. The port strips the prefix first, so the rule does fire. Expect differences on multi-cluster tumours if the old run used prefixed labels.
- **G: a plain cluster label doesn't disallow its own clone** (`4` is listed, not `clone4`). A ploidy-only event can therefore land on the same clone as a mutation-supported one. Kept as is.
- **H and K1** would raise on `alpaca_None`, on unresolved multi-cluster labels (`3,6`), or when a cluster is missing from `mut_counts`. The port skips these, or gives a null count.
- **G stale variable** (section 5, item 4): when there are no other events, the disallowed ancestors are empty.
- **class / sub_class:** class maps the triage (no_gd → no_wgd, clonal_gd_only → clonal_wgd, alpaca_called → ploidy_only, mut_supported). sub_class isn't in the notebook. It follows the tx_data `wgd_calls` definitions and is computed only when every event is placed on an ALPACA clone.
- **ALPACA score** defaults to `ratio_A` (notebook). `--alpaca-score max_AB` applies decision 2.

## 5. Things probably not to replicate

1. **Greedy, one-event-at-a-time ALPACA assignment (C).** Events take the best remaining candidate clone in file order. Nothing requires the chosen clone to explain *all* regions of the event and *none* of the undoubled regions, and nothing orders a first GD above a second GD.

   *Revised 2026-11-09 after Piotr's review.* The first version of this example was wrong on two counts:
   - Its ratios were inconsistent. ALPACA ratios are relative to the **parent** (`calculate_wgd_ratios`), so c2 = 1.95 and c4 = 1.6 meant two sequential doublings on the c2→c4 lineage, and nothing explained c5.
   - PGDD would not produce two events. With no doubled clusters, it puts all regions with the same ploidy GD status into one event (`seperate_gd_events`, `R/parrallel_gd_caller.R` ~l.490). Checked by running `detect_par_gd` on that input: one event, `samples = R1,R2`.

   Piotr's parsimony reading is the right answer for that tree: **one GD on c2**. It can't be on the MRCA because c3 sits in an undoubled region (R3), and c2 explains both c4 and c5. Tree + region ploidy alone pin it down, without ratios.

   Where the old rule actually fails is with nested (first + second) ploidy-only events:

   ```
   c1 (MRCA) ── c2 ──┬── c4      R1 = 100% c4   num_gds 2
            │        └── c5      R2 = 100% c5   num_gds 1
            └── c3               R3 = 100% c3   num_gds 0
   ALPACA max(A,B) vs parent: c2 = 1.95, c4 = 2.05, c5 = 1.0, c3 = 1.0
   PGDD: E1 (first GD) in R1,R2; E2 (second GD) in R1
   ```
   The truth is GD1 on c2 and GD2 on c4. The old rule picks E1 from {c2, c4, c5} → **c4** (2.05 > 1.95), then E2 from {c4, c2} minus c4 → c2. That gives the first GD on c4, which doesn't cover R2, and the second GD on its ancestor. Both placements are impossible, yet the tumour is called `resolved` (both > 1.5). Processing E2 first happens to give the right answer, so order matters as well.

   **Status:** deferred at first, so the pipeline could be tested on real data with the old rule (C). **Replaced on 2026-11-09** by the tumour-level ploidy fill, after the real cases LTX0025 and LTX0073 (Piotr): see "Ploidy fill" below. It covers both C and G.

   Concerns with the current rule, in short:
   - **Coverage is not checked.** The chosen clone may not be ancestral to all of the event's regions.
   - **Exclusion is not checked.** The clone may be ancestral to an undoubled region; only the MRCA is excluded explicitly.
   - **GD order is not enforced.** A second GD can land above the first.
   - **Results depend on order.** Events take the best remaining clone in file order.
   - **One event, one clone.** Events are never merged or split, and the clone is picked by ratio even when tree + ploidy already determine it.
   - **The resolved check looks at ratios only.** A tumour can be `resolved` (all ratios > 1.5) with impossible placements, as in the example above.

   **Proposed rule (later improvement):** first compute each event's *feasible* clones from tree + region ploidy alone:
   - the clone must be an ancestor-or-self of every clone in the event's regions;
   - it must not be an ancestor-or-self of any clone in a region without the event;
   - a second GD must be a descendant-or-self of the first GD's clone on that lineage.

   Ratios are used only to choose within a feasible set that has more than one clone (a chain of single-child clones). If the feasible set is empty, flag the event rather than force it.
2. **The MRCA is excluded from ALPACA candidates (C).** A ploidy-only "subclonal" event can never land on the MRCA, even when the MRCA has the highest `ratio_A`.
3. **The multi-cluster rule contradicts its comment (E).** The comment says parent; the code picks the child.
4. **Stale variable (G).** `forbidden_descendants` is reset for each event, but `forbidden_ancestors` is only set when other events exist. If a tumour's ploidy-only event is its only event, the code either raises a `NameError` (if no earlier event set the variable) or silently reuses the **previous tumour's** disallowed ancestors. *Latent in the old run:* G only sees events of mutation-supported tumours, which always have another event. It would bite as soon as G is applied more widely.
5. **The ancestor check doesn't do what its comment says (G).** The comment says "check if all its children are present in the sample", but the code only checks for disallowed descendants. An ancestor can be chosen even when its descendants sit in samples that are **not** doubled. That contradicts the region pattern of the event.
6. **The ploidy check compares only maxima (H).** Per-sample consistency isn't checked, and it relies on the cohort summary when `GDs_per_region.num_gds` would do.
7. **The reversion check (I) uses allele A only, by median, and only in mutation-supported tumours.**
8. **The quality rule (K)** uses a cohort-relative 5th-percentile cut-off, so it changes with the cohort and the run. It also makes a whole tumour Low because of one event.
9. **State tracking through string prefixes** (`alpaca_`, `clonal_`, `multicluster_`, `ploidy_`, `false_call`, `reversion_`) packed into `clusters`, plus several `events_df` copies. We'll use explicit columns (`placed_clone`, `placement_method`, `flags`).
10. **Hard-coded tumour lists, side-effect writes to shared queue files, and `which_q` lookups.** These are interactive scaffolding.
11. **`gd_p_value` and the hypergeometric 2×2 tables.** We don't use the hypergeometric test any more. The current run keeps only bootstrap-stable calls (section D).

### Ploidy fill (step 04 default since 2026-11-09; `--ploidy-rule notebook` restores rules C and G)

This replaces the per-event rules C and G with one tumour-level rule. PGDD's ploidy-only events are not placed one by one. They are **superseded**: GDs are added until every sample's GD level matches its PGDD `num_gds`.
1. **Fixed GDs:** clonal events go on the MRCA, and mutation-supported events go on their clone (multi-cluster parent–child pairs go to the child). If a fixed event can't be put on an ALPACA clone (off-tree cluster, unresolved multi-cluster), the tumour is left unresolved (`fill = fixed_event_not_on_tree`).
2. **Sample level:** a clone's count is the number of GDs on its root path. A sample's level is the largest k such that clones with count ≥ k make up **more than 50%** of the sample, from the ALPACA `cp_table` (the sample-level analogue of PGDD's ≥ 50%-of-genome rule). Targets are `num_gds` from `pgdd_GDs_per_region`.
3. **Adding GDs:** one at a time, on the clone that (a) raises the most samples still below target, without pushing any sample above its target; (b) leaves the fewest cells **above** their sample's target (no over-doubling of clones in samples that already have enough); (c) leaves the most cells **at** their sample's target (prefers the ancestral clone that also doubles minor clones of already-doubled samples); (d) is the **most ancestral** (fewest ancestors), so when the data can't tell clones apart the GD goes as high as allowed (added 2026-11-10, LTX1338); (e) has the highest ALPACA score; (f) comes first in the tree. `fill_tie` now means only the ALPACA score separated equally deep clones. Cell fractions are summed over samples. (Revised 2026-11-10 after LTX0501. The first version used "fewest cells in samples not raised", which picked clone9 there.) A clone can take several GDs. Each added GD is a row `fill<i>` in `events_placed.tsv`, with the samples it raised and any equally good alternatives (`fill_alternatives`, flag `fill_tie`).
3b. **Relaxed pass** (2026-11-11, LTX1055): if a sample is still short and no clone is the majority of it without over-doubling other samples, the clone proportions and the sample's ploidy disagree. The GD then goes on the clone covering the largest share (≥ 10%) of the still-short samples, as long as no cell anywhere ends above its sample's target (tie-breaks: ancestral, then ALPACA score). It's flagged `minority_coverage(<share>:<sample>)` and penalised heavily: factor 0.2 at 0% coverage, rising to 0.5 at ≥ 50% (Piotr, 2026-11-11: place it for transparency, but penalise heavily).
4. **Status:** `resolved` when every sample ends exactly at its `num_gds` (`fill = complete`). Otherwise `unresolved`, with `fill = missing:<sample>=n` or `over:<sample>=n` (the fixed GDs already exceed the target).
5. The H check (GDs on a path vs sample ploidy) counts events under this rule.

**LTX0025-T1:** clonal GD on clone1; R4 needs two more. R4 = clone3 (0.6) + clone5 (0.4). clone2 would also raise R2 (its subtree is all of R2) past its single GD. clone5 alone is under 50% of R4. So clone3 takes both GDs → clone1 = 1, clone3 = 2, resolved. The old rule gave clone3 then clone2 (unresolved, fixed by hand).

**LTX0073-T1** (mutation + ploidy): fixed GDs on clone1 (clonal) and clone3 (mutation-supported, cluster 3). Levels then: R1 = 1 (target 3), R2 = 2 (2), R3 = 2 (3), R4 = 2 (2), because clone3's subtree is 63% of R4. PGDD's ploidy-only event GD3 (R1 + R4) can't fit this tree, which is why per-event placement fails. Fill: clone7 twice (R1 = clone7 + clone15), then clone14 (R3 = clone14 + clone16). clone5 would also raise R3, but it covers clone12, 30% of R4: more spill. Result clone1 = 1, clone3 = 1, clone7 = 2, clone14 = 1, resolved. That's Piotr's answer.

**LTX0501-T1** (2026-11-10): mutation-supported GD on clone4 (R1, 92% clone4 + clone3); PGDD's ploidy-only GD1 covers R2 + R3, with no clonal GD. No single clone fits R2 + R3: clone1 would also raise R1 past 1. Fill: R3 → clone5 (the only option; clone2 is 39%). R2 → clone6 or clone9 both raise it with no over-doubling, but clone6 also brings its 8% share of R1 to R1's target (1 GD). Under (c) that's clone6. Result clone4, clone5, clone6, resolved: the old hand-curated answer.

**LTX1338-T1** (2026-11-10): mutation-supported GD on clone2 (R1, R2). It also explains R3 through clone11, 65% of R3, which is why PGDD's ploidy-only GD2 (R3 + R4) is only needed for R4. R4 = clone4 (24%, under clone2) + clone15/clone14/clone5 (76%, under clone9). clone9 (extinct) and its only child clone15 cover exactly the same cells, and ALPACA can't separate them (clone15's ratio is exactly 1.0). The old order fell through to the ALPACA score (ratio_A clone9 0.91 < clone15 1.0) and picked clone15. With (d), the GD goes on clone9, Piotr's answer.

**LTX1055-T1** (2026-11-11): mutation-supported GD on clone4 (brain met FR1, FR2). Primary R1 is fully doubled (`genome_frac_nMaj2` 0.9997, ploidy 3.9), but `cp_table` gives it 68% clone2, the same clone that makes up all of the undoubled R2–R7, and only 32% clone3. Doubling clone2 (or clone1) would push R2–R7 above 0, and clone3 alone is short of 50% of R1, so the strict fill stopped with R1 missing a GD (unresolved). The relaxed pass puts it on clone3, which ALPACA supports (ratio 3.2 vs clone2). Result: clone3 + clone4, resolved, with a heavy `minority_coverage(0.32:R1)` penalty (quality ≈ 0.39). The R1 clone proportions are probably underestimating clone3.

**LTX1018-T1** (2026-11-11): clonal GD on clone1 and mutation-supported GD on clone3 (PGDD cluster 3). Primary R2 needs a third GD. clone2 (+ clone7) is only 34% of R2. clone5 is 21% of R2 but would also double clone8 and its descendants in R1, which has 1 GD (over-doubling, so not allowed; the reversion reading isn't plausible either). No clone explains R2 cleanly, so the relaxed pass forces the GD onto clone2 with a heavy penalty: event quality ≈ 0.3 (`minority_coverage(0.34)`, plus borderline and ambiguity factors).

All six cases were replayed together from the real PGDD files.

## 5b. ALPACA ratio threshold: concerns (Darwin, 2026-11-09; revisit later, 1.5 in use)

`ratio_X` = the length-weighted mean of child CN / parent CN over segments where the parent allele is > 0. The MRCA is compared with a diploid 1|1. A clean doubling gives 2 and no change gives 1, so 1.5 is the natural midpoint. **Decision: use 1.5 for now** (Piotr).

Concerns:
1. **High-level gains mimic a doubling.** The mean of ratios is not robust. A non-doubled allele unchanged over 90% of the genome and at 6× the parent over 10% scores 0.9 + 0.6 = 1.5. Arm gains or amplifications over ~10–20% of the genome can cross the threshold without a WGD.
2. **Losses after a doubling hide it.** A doubled allele that falls back to one copy over half the genome scores 1.5. `max(A,B)` mitigates this when the other allele kept its copies.
3. **`max(A,B)` is biased upward when one allele has little data.** If the parent has LOH over much of the genome, that allele's ratio rests on few segments (low `weight_sum_X`), and a single gained segment can push it over 1.5. Idea: let an allele into the max only if its `weight_sum` covers enough of the genome (e.g. ≥ 30%), otherwise use the other allele or `ratio_total`.
4. **A hard cut.** Values near the threshold (say 1.3–1.7) are weak evidence. Idea: keep 1.5 for the call, but add `alpaca_ratio_borderline` to the quality.
5. **Inherited errors.** The ratio is relative to the parent's ALPACA CN, so an error in the parent's CN propagates to all its children.

Checks once real data exist:
- The distribution of `max(A,B)` across clones should be bimodal near 1 and 2, with the trough near 1.5.
- Clones carrying mutation-supported PGDD GDs should mostly be > 1.5; clearly undoubled clones (diploid regions) mostly < 1.5.

## 6. Decisions (answered by Piotr 2026-11-07 unless marked open)

1. **QC (A).** ✅ Keep both the manual-QC `FINALIZED` filter and the truncal-CCF > 90% filter, sourced from tx_data.
2. **ALPACA `wgd_ratio_scores.csv`.** ✅ Generated by ALPACA (`alpaca wgd`); in tx_data as `alpaca_wgd_ratios` (draft PR #7). Path `tracerx_842/ALPACA/output/tumour_results/<tumour>/wgd_ratio_scores.csv`, 12 columns, consistent with `all_tumours_combined.csv` (Piotr, 2026-11-09). **Score = `max(ratio_A, ratio_B)`** (2026-11-09): after a doubling, tumours accumulate losses, which show in the minor allele; with phased data B is not always the minor allele. The 1.5 threshold for the new ALPACA is still open.
3. **Multi-cluster events (E).** ✅ Parent–child pair → the **child**, with a quality penalty (Piotr, 2026-11-09; reasoning in E). Non-adjacent pairs and > 2 clusters: open.
4. **Ploidy-only placement (C vs G, multi-sample events).** ⏸ Later, decided on a few worked examples (reminder 3).
5. **wgd_viewer in the loop / step 06 format.** ⏸ Last.
6. **Second GDs.** Open (they were skipped entirely in G).
7. **Triage classes (B).** Open: keep the four, or the placeholder's `clonal_wgd / mut_supported / ploidy_only / no_wgd` plus our `ploidy_only_pruned`.
8. **Quality (K).** Open: draft proposal in K'.

## 7. Open checks and reminders

| # | What | Who | Status |
|---|---|---|---|
| 1 | Multi-cluster events (E): parent or child? | Piotr | ✅ child + quality penalty (2026-11-09) |
| 2 | ALPACA ratio threshold: is 1.5 still right for the new ALPACA, now applied to `max(ratio_A, ratio_B)`? | Piotr | score ✅ max(A,B); threshold open |
| 3 | Ploidy-only placement rule: pick a few real tumours and decide C vs G, and multi-sample events. | Piotr + Darwin | later |
| 4 | Quality thresholds `N_MIN`, `M_MIN`. (Source of the per-sample doubled-genome fraction ✅: steps 00/03.) | Piotr | open |
| 5 | Greedy ALPACA assignment (section 5, item 1). | Piotr + Darwin | ✅ replaced by the tumour-level ploidy fill, 2026-11-09 (LTX0025, LTX0073) |
| 6 | ALPACA threshold concerns (section 5b). | Piotr | revisit later; 1.5 in use |
