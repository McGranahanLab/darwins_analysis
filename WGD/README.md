# WGD

TRACERx whole-genome doubling: build the **integrated WGD table** by
combining ParallelGDDetect (PGDD) WGD calls with ALPACA clone copy number
and the CONIPHER tree. Requested by Piotr (2026-11-06). WES cohort first,
WGS later.

The `wgd_calls` table currently in tx_data is a made-up placeholder. This
analysis replaces it.

## Inputs

| Source | Where | Notes |
|---|---|---|
| PGDD outputs | `/nemo/.../tracerx_842/WGD/release_tx842/pgdd` (mock: tx_data PR #5) | **Old run.** Will be regenerated with the new hg19 ALPACA results, then again for hg38. Only stability-filtered `pgdd_*` files are read; `unstable_calls_*` is for manual review only. |
| ALPACA clone CN | tx_data `alpaca` | |
| Tree clones | tx_data `clone_proportions` (ALPACA `cp_table.csv`) | The tree is the set of clones ALPACA has. PGDD runs on pre-tree clusters, so some called clusters are not on it. |

PGDD: amf71/ParallelGDDetect master `cb0437b` (bootstrap). IDs are
canonical: `tumour_id` = `LTX0001-Tumour1`, `sample_id` = muttable
`sample_name_hash`, `cluster_id` = PGDD cluster as a string (`"3"`, which is
`clone3` on the tree).

## Steps

| # | Step | Status |
|---|---|---|
| 01 | Parse PGDD; explode comma-joined `clusters`, `samples`, `gd_clusters` into long tables | ✅ `scripts/01_parse_pgdd.py` |
| 02 | Prune against the tree (= clones ALPACA has: `cp_table.csv` / `tree_paths.json`, via tx_data `clone_proportions`): drop tumours without a tree; remove off-tree supporting clusters; an event left with no cluster is kept as `ploidy_only_pruned` | ✅ `scripts/02_prune_tree.py` |
| 03 | Place events on the clone tree (event, regions and clusters to a clone/branch) | ⏳ |
| 04 | Check against ALPACA clone CN | ⏳ |
| 05 | Integrated per-clone WGD table | ⏳ |

## Run

```bash
pixi install
pixi run parse_pgdd --nemo-root /nemo --out output/<name>      # on Nemo
pixi run parse_pgdd --nemo-root ../../tx_data/nemo_mock --out output/mock
pixi run prune_tree --tx-data-dir ../../tx_data/data --out output/mock
pixi run test     # needs tx_data with PR #5 next to darwins_analysis, or TX_DATA_REPO
```

Step 01 output (`output/<name>/01_pgdd/`):

- `pgdd_event_clusters.parquet`: one row per event × supporting cluster.
  Events without mutation support (clonal or ploidy-only) keep a single row
  with a null `cluster_id`.
- `pgdd_event_regions.parquet`: one row per event × region.
- `pgdd_region_clusters.parquet`: one row per region × doubled cluster.
- The raw tables as parquet, plus `run.json` (bootstrap settings and counts).

Step 02 output (`output/<name>/02_pruned/`):

- `events_pruned.parquet`: one row per kept event. It has the PGDD columns plus
  `clusters_pruned`, `clusters_removed`, `n_clusters_raw`, `n_clusters_kept`
  and `support`. `support` is one of `clonal`, `mut_supported`, `ploidy_only`
  or `ploidy_only_pruned`. `ploidy_only_pruned` means the regions carry the
  event, but its clone is unknown because all supporting clusters were off
  the tree. PGDD's own `is_subclonal_mutation_supported` is left as written.
- `events_dropped.parquet`: events in tumours with no tree (`drop_reason`).
- `event_clusters_flagged.parquet`: event × cluster rows with `clone` and `on_tree`.

Matching regions across sources uses the sample stem (`LTX0001_SU_T1-R1`),
not the hash suffix.

## Docs

- [`docs/previous_analysis_logic.md`](docs/previous_analysis_logic.md): the logic of the previous consolidation (`templates/previous_analysis.ipynb`) as a pipeline, with problems and open decisions.

## Caveats

- PGDD groups clusters into events with `grepl(cluster, gd_clusters)`
  (`R/parrallel_gd_caller.R:412`), which is a substring match: cluster `2`
  also matches `12`. In the mock, LTX0007's parallel GDs on clusters 2 and
  12 are merged into one event. Piotr is looking at this upstream.
