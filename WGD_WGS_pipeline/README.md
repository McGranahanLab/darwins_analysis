# WGD_WGS_pipeline

PGDD → WGD caller → ALPACA post-hoc correction on the TRACERx **WGS** cohort, as one bash pipeline with one
config file (Piotr, 2026-11-16). Inputs come from the tx_data store; the WGD caller and the correction are
the code in [`../WGD`](../WGD/) and [`../alpaca_post_hoc_correction`](../alpaca_post_hoc_correction/),
called as they are.

```bash
./run.sh config.toml                  # everything
./run.sh config.toml --from wgd04     # resume (steps: check stage wgd00 pgdd wgd01 wgd02 wgd03 wgd04 wgd05 wgd07 correct)
./run.sh config.toml --only pgdd
```

## Steps

| step | what | writes (under `run.out_dir`) |
|---|---|---|
| check | config valid; every tx_data table present (lists the missing ones and stops) | `logs/check.log` |
| stage | tx_data WGS tables → the layouts the steps read; `[pgdd.settings]` → JSON, `[wgd.params]` → TOML | `inputs/`, `config/` |
| wgd00 | WGD step 00: per-sample `num_gds`, PGDD input; single-clone tumours set aside | `wgd/00_pgdd_input/` |
| pgdd | `detect_par_gd()` with `[pgdd.settings]` (`pipeline/run_pgdd.R`) | `pgdd/` |
| wgd01–03 | parse PGDD, prune to the ALPACA tree, per-sample summary | `wgd/01_pgdd/` … `03_summary/` |
| wgd04 | place GDs with `[wgd.params]`; QC if `[qc].enabled` | `wgd/04_placed/` |
| wgd05 | per-clone table; GDs below `[wgd].quality_cutoff` rejected | `wgd/05_table/` |
| wgd07 | cohort report | `wgd/07_report/report.html` |
| correct | ALPACA post-hoc correction on the filtered calls | `correction/` |

Every run copies the config to `logs/config.used.toml` and records the darwins_analysis commit; each step
logs to `logs/<step>.log`.

## Config (`config.toml`)

One file for all steps: `[run]` (output dir, tumour subset), `[tx_data]` (store), `[pgdd]` (checkout,
R command) and `[pgdd.settings]` (passed verbatim to `detect_par_gd`), `[wgd]` (**`quality_cutoff`**,
mutation-supported rejection, ALPACA score, report threshold), `[qc]`, `[correction]`, and `[wgd.params]`,
the WGD caller's rules (same keys and defaults as `WGD/params.toml`; a test keeps them equal).

## Inputs from tx_data

| table | used by | source on Nemo |
|---|---|---|
| `wgs_consensus_segments` | step 00 (num_gds) | Battenberg `tree_input/` |
| `wgs_clone_proportions` | steps 02, 04, correction | `wgs_resources/alpaca/input/<label>/cp_table.csv` |
| `wgs_conipher_mutations` | step 00 (SCoutput.CLEAN + treeTable) | CONIPHER `2026_08_mean_max_absent/.../snv/<label>/` |
| `wgs_conipher_tree_edges` | step 00 (default tree, `allTrees.txt`) | same |
| `wgs_conipher_cluster_info` | step 04 CCF QC (`clusterInfo.txt`) | same |
| `wgs_alpaca` | step 04, correction (`ALPACA_output_<label>.parquet`) | `wgs_resources/alpaca/output/tumour_results/<label>/` |
| `wgs_alpaca_wgd_ratios` | steps 02, 04, 05, correction (`wgd_ratio_scores.parquet`) | same |

The last five came with tx_data PR #11 (2026-11-17). **Manual QC** is not a tx_data table: the WGS
run reads the WES manual QC file (`[qc].manual_qc_table`; Piotr, 2026-11-17: every WGS tumour is in
the WES cohort). tx_data's `nemo_mock/` has a minimal mock of it at the mirrored path.

## Sandbox

```bash
cd ../../tx_data && pixi run build sandbox && cd -   # tx_data mock cohort, with tx_data #11
./run.sh config.sandbox.toml                          # writes output/sandbox/
```

On the mock (2026-11-17, with the single-clone change): all 11 tumours reach `wgd_calls.tsv`. The
8 single-clone ones are placed by ploidy: LTX0001 and LTX0102 (WGD-like, ploidy 3.6–3.9) get one
clonal GD, the rest none, all quality 1.0 (no borderline sample). LTX0004 gets a clonal GD from PGDD;
LTX0103 a borderline clonal GD (genome fraction 0.53, quality 0.07) that the 0.3 cut-off rejects;
LTX0104-Tumour2 is `qc_failed` (truncal CCF 86) and LTX0105 too (not FINALIZED in the mock QC
table). The correction runs on all 11 (none changed). No PGDD settings warnings.

## Requirements

`pixi` (runs the WGD and correction environments), R with `data.table` + `jsonlite` (`[pgdd].rscript`,
e.g. `micromamba run -n pgdd Rscript`), and an amf71/ParallelGDDetect checkout at `cb0437b` (`[pgdd].src`).

## Tests

```bash
pixi run --manifest-path ../WGD/pixi.toml pytest tests/ -q
```

`tests/make_fixture.py` writes a synthetic WGS store with every table (including the proposed ones);
the last test runs the whole pipeline on it, PGDD included (skipped without pixi + micromamba).

## Single-clone tumours

Step 00 sets single-clone tumours aside (they don't go to PGDD). Since 2026-11-17 (Piotr) step 04 reads
them back (`00_pgdd_input/single_clone_samples.tsv`): triage `single_clone`, no events, and the ploidy
fill places GDs on their only clone until every sample matches its `num_gds`, with the ploidy prior,
quality factors and QC as for any tumour. Class `clonal_wgd` or `no_wgd`. Samples that disagree and
can't be reconciled leave the tumour `unresolved`.

## Known gaps / TODO

- **TODO: a separate WES pipeline** (Piotr, 2026-11-19). This pipeline is strictly WGS: its PGDD thresholds
  are WGS-scale (`discover_num_muts_threshold` 100, `discover_num_2_cpn_muts_threshold` 50,
  `check_num_2_cpn_muts_threshold` 30). Create `WGD_WES_pipeline/`, identical in structure, with its own
  config (WES settings, WES tx_data tables), so WES settings are never used for WGS or the reverse.
- TODO (tx_data): give the manual QC mock the real file's full column set.
- `check` stops when `[pgdd].src` is a PGDD without the bootstrap/stability arguments (e.g. the wlippa
  fork at 1af16a8); use amf71/ParallelGDDetect at cb0437b or later.

## Output: quality_pass

`wgd_calls.tsv` keeps every call. `quality_pass` (2026-11-19) is False on a clone carrying a GD with event
quality below `[wgd].quality_cutoff` (the GDs `wgd_calls_rejected_below_<Q>.tsv` drops; mutation-supported GDs
pass unless `reject_mutation_supported`), True otherwise, null on rows without calls or with cutoff -1.
