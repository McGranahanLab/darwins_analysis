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
2. **Truncal CCF QC.** For each passing tumour, read CONIPHER `clusterInfo.txt` and keep clusters with `treeClust & truncal`. Take the **lowest** `meanCCF` among truncal clusters. The tumour passes if that value is above **0.9**.
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
3. Pick the candidate with the **highest `ratio_A`**. Write `clusters = "alpaca_<clone>"` and store `alpaca_ratio_A`.

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
  - **Decision (Piotr, 2026-11-07): follow the code (child) for now.** ⏰ Reminder 1 (section 7): Piotr checks real output to establish which version is correct.
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

#### K'. Proposed replacement quality (draft, 2026-11-07)
Piotr: the metric needs improving, with three components. Proposal: compute each component **per event**, keep them as columns, and derive the label from fixed thresholds (not cohort percentiles).

| Component | Per-event value | Flag when | Source |
|---|---|---|---|
| **Mutation support** | `n_support` = `num_cn2` of the supporting cluster (max over the event's regions); also `stability` / `q_value` from the bootstrap summary | `n_support < N_MIN` (e.g. 10, TBD). Ploidy-only events have no mutation evidence and are flagged `no_mut_support` | `pgdd_mut_counts`, `pgdd_bootstrap_discovery_summary` |
| **Borderline ploidy** | `ploidy_margin` = min over the event's regions of (fraction of genome with major CN ≥ 2) − 0.5. For a second GD the fraction uses major CN ≥ 3. This is the "Bielski fraction" of `[48]` | `ploidy_margin < M_MIN` (e.g. 0.05–0.10, so 55% doubled is borderline). Arguably also check the regions *without* the event: 45% is a borderline "no". | Not in the PGDD outputs: PGDD takes `num_gds` as an input (threshold ≥ 50% with major CN ≥ 2, `R/parrallel_gd_caller.R:45–50`). Needs the per-sample fraction from the PGDD input preparation or the ASCAT segments. |
| **GD order** | `gd_order` = 1 for a first GD, 2 for a second, 3 for a third, counted along the root→leaf lineage | `gd_order ≥ 2`, stronger penalty for 3 | Our placement (03) |

Combining: event quality = `High` with no flags, `Medium` with one, `Low` with two or more (or with `no_mut_support` plus `borderline_ploidy`). Tumour quality = the worst event, but the per-event flags are kept so users can filter on the reason. Thresholds `N_MIN` and `M_MIN` are open (section 7, reminder 4).

## 4. How it maps to the current WGD steps

| Previous | Current |
|---|---|
| A. QC | not yet. It needs the QC table and CONIPHER `clusterInfo`, or tx_data equivalents. |
| (implicit) | **02** pruning: off-tree clusters and tumours without a tree. The notebook never handles off-tree clusters explicitly. `get_descendants` on an off-tree clone just returns `[]`. |
| B. triage | 03, derived from `events_pruned.support` rather than PGDD's per-tumour flags |
| E, F, G, C | 03 placement |
| H, I | 04 checks against ALPACA |
| J, K | 05 integrated table |

## 5. Things probably not to replicate

1. **Greedy, order-dependent ALPACA assignment (C).** Events take the best remaining clone in file order, so a different order can produce a different assignment. Synthetic example:

   ```
   c1 (MRCA) ── c2 ──┬── c4      R1 = 100% c4
            │        └── c5      R2 = 100% c5
            └── c3               R3 = 100% c3   (not doubled)
   ALPACA ratio_A: c2 = 1.95, c4 = 1.6, c5 = 1.1, c3 = 1.0
   PGDD: two ploidy-only events, E1 in R1 and E2 in R2 (1 GD per region)
   ```
   Candidates: E1 → {c4, c2}; E2 → {c5, c2} (sample clones plus ancestors, MRCA excluded). c2 is the best candidate for both.
   - Order E1, E2: E1 → c2 (1.95), E2 → c5 (1.1). Tumour **unresolved** (1.1 < 1.5).
   - Order E2, E1: E2 → c2 (1.95), E1 → c4 (1.6). Tumour **resolved**.

   So the same tumour flips between resolved and unresolved with row order. Both answers are also biologically inconsistent: one event on c2 and the other on its descendant means one region carries **two** GDs, while PGDD saw one per region. Here ALPACA explained both doublings parsimoniously with one GD on c2. The two consistent answers are "one shared GD on c2" (PGDD's two events are really one) or "parallel GDs on c4 and c5". A greedy one-event-at-a-time pick can produce neither. A joint assignment (all events at once, no nesting allowed between first GDs, maximise total score) would.
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

## 6. Decisions (answered by Piotr 2026-11-07 unless marked open)

1. **QC (A).** ✅ Keep both the manual-QC `FINALIZED` filter and the truncal-CCF > 0.9 filter, sourced from tx_data.
2. **ALPACA `wgd_ratio_scores.csv`.** ✅ Generated by ALPACA (`alpaca wgd`); add it to tx_data. In progress (tx_data draft PR). The 1.5 threshold for the new ALPACA is still open.
3. **Multi-cluster events (E).** ✅ For now follow the code (child), not the comment. Reminder 1. Non-adjacent pairs and > 2 clusters: open.
4. **Ploidy-only placement (C vs G, multi-sample events).** ⏸ Later, decided on a few worked examples (reminder 3).
5. **wgd_viewer in the loop / step 05 format.** ⏸ Last.
6. **Second GDs.** Open (they were skipped entirely in G).
7. **Triage classes (B).** Open: keep the four, or the placeholder's `clonal_wgd / mut_supported / ploidy_only / no_wgd` plus our `ploidy_only_pruned`.
8. **Quality (K).** Open: draft proposal in K'.

## 7. Open checks and reminders

| # | What | Who | Status |
|---|---|---|---|
| 1 | Multi-cluster events (E): check real output to establish whether the event belongs on the parent or the child. Until then the code's behaviour (child) is used. | Piotr | ⏰ open |
| 2 | ALPACA ratio threshold: is 1.5 still right for the new ALPACA (and now that `ratio_B` / `ratio_total` exist)? | Piotr | open |
| 3 | Ploidy-only placement rule: pick a few real tumours and decide C vs G, and multi-sample events. | Piotr + Darwin | later |
| 4 | Quality thresholds `N_MIN`, `M_MIN`, and a source for the per-sample doubled-genome fraction. | Piotr | open |
