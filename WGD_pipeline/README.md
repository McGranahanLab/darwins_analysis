# WGD_pipeline

PGDD → WGD caller → ALPACA post-hoc correction on TRACERx, as one bash pipeline (Piotr, 2026-11-16) with
one config file per cohort (2026-11-24; until then `WGD_WGS_pipeline/`, WGS only):

| config | cohort | PGDD thresholds (discover_num_muts / discover_num_2_cpn / check_num_2_cpn) | tx_data tables |
|---|---|---|---|
| `config.wgs.toml` | WGS | 100 / 50 / 30 (Piotr, 2026-11-19) | `wgs_*` |
| `config.wes.toml` | WES (tx842) | 10 / 5 / 3 (Piotr, 2026-11-09; = `wgd.pgdd.EXPECTED_SETTINGS`) | unprefixed (`alpaca`, …) |
| `config.wgs.sandbox.toml` | WGS, tx_data mock | as WGS | `wgs_*` |

Inputs come from the tx_data store; the WGD caller and the correction are the code in [`../WGD`](../WGD/)
and [`../alpaca_post_hoc_correction`](../alpaca_post_hoc_correction/), called as they are. The steps don't
know the cohort: `[tx_data.tables]` says which tx_data table feeds each input, and `stage` writes the same
layout for both.

```bash
./run.sh config.wgs.toml                  # everything, WGS (writes output/wgs/)
./run.sh config.wes.toml                  # everything, WES (writes output/wes/)
./run.sh config.wgs.toml --from wgd04     # resume (steps: check stage wgd00 pgdd wgd01 wgd02 wgd03 wgd04 wgd05 wgd07 correct)
./run.sh config.wgs.toml --only pgdd
```

**Keeping WES and WGS apart.** `check` stops if `[run].cohort` and `[tx_data.tables]` disagree (a WGS config
naming a non-`wgs_` table, or the reverse). A test (`test_wes_and_wgs_configs_differ_only_in_cohort_settings`)
keeps the two configs equal except for `[run].cohort`, `[run].out_dir`, `[tx_data.tables]` and the three PGDD
mutation-count thresholds, so a change to the shared settings (`[wgd]`, `[wgd.params]`, …) must be made in both.

## Steps

| step | what | writes (under `run.out_dir`) |
|---|---|---|
| check | config valid; every tx_data table present (lists the missing ones and stops) | `logs/check.log` |
| stage | the cohort's tx_data tables → the layouts the steps read; `[pgdd.settings]` → JSON, `[wgd.params]` → TOML | `inputs/`, `config/` |
| wgd00 | WGD step 00: per-sample `num_gds`, PGDD input; single-clone tumours set aside | `wgd/00_pgdd_input/` |
| pgdd | `detect_par_gd()` with `[pgdd.settings]` (`pipeline/run_pgdd.R`) | `pgdd/` |
| wgd01–03 | parse PGDD, prune to the ALPACA tree, per-sample summary | `wgd/01_pgdd/` … `03_summary/` |
| wgd04 | place GDs with `[wgd.params]`; QC if `[qc].enabled` | `wgd/04_placed/` |
| wgd05 | per-clone table; GDs below `[wgd].quality_cutoff` rejected | `wgd/05_table/` |
| wgd07 | cohort report | `wgd/07_report/report.html` |
| correct | ALPACA post-hoc correction on the filtered calls | `correction/` |

Every run copies the config to `logs/config.used.toml` and records the darwins_analysis commit; each step
logs to `logs/<step>.log`.

## Config (`config.wgs.toml`, `config.wes.toml`)

One file for all steps: `[run]` (cohort, output dir, tumour subset), `[tx_data]` (store) and `[tx_data.tables]`
(pipeline input → tx_data table), `[pgdd]` (checkout,
R command) and `[pgdd.settings]` (passed verbatim to `detect_par_gd`), `[wgd]` (**`quality_cutoff`**,
mutation-supported rejection, ALPACA score, report threshold), `[qc]`, `[correction]`, and `[wgd.params]`,
the WGD caller's rules (same keys and defaults as `WGD/params.toml`; a test keeps them equal).

## Inputs from tx_data

| input (`[tx_data.tables]` key) | used by | WGS table (source on Nemo) | WES table |
|---|---|---|---|
| `segments` | step 00 (num_gds) | `wgs_consensus_segments` (Battenberg `tree_input/`) | `msai_segments` — **not in tx_data yet** |
| `clone_proportions` | steps 02, 04, correction | `wgs_clone_proportions` (`wgs_resources/alpaca/input/<label>/cp_table.csv`) | `clone_proportions` |
| `conipher_mutations` | step 00 (SCoutput.CLEAN + treeTable) | `wgs_conipher_mutations` (CONIPHER `2026_08_mean_max_absent/.../snv/<label>/`) | `conipher_mutations` — **not in tx_data yet** |
| `conipher_tree_edges` | step 00 (default tree, `allTrees.txt`) | `wgs_conipher_tree_edges` (same) | `conipher_tree_edges` — **not in tx_data yet** |
| `conipher_cluster_info` | step 04 CCF QC (`clusterInfo.txt`) | `wgs_conipher_cluster_info` (same) | `conipher_cluster_info` — **not in tx_data yet** |
| `alpaca` | step 04, correction | `wgs_alpaca` (`wgs_resources/alpaca/output/tumour_results/<label>/`) | `alpaca` |
| `alpaca_wgd_ratios` | steps 02, 04, 05, correction | `wgs_alpaca_wgd_ratios` (same) | `alpaca_wgd_ratios` |

**Why `segments` comes from different places (Piotr, 2026-11-24).** Both cohorts need fractional allele-specific
CN per segment per sample, MSAI-corrected, from the same solution that feeds CONIPHER and ALPACA. WES runs ASCAT,
then Refphase; Refphase detects MSAI itself, so its output is already MSAI-corrected and there are no separate
`MSAI_adjusted` files. The WES `segments` therefore come from the ALPACA input (same Refphase solution as CONIPHER).
WGS runs Battenberg, which does no MSAI detection, so a post-hoc MSAI correction produces the
`*_MSAI_adjusted_*segments_tree_input.tsv` files. The table name `msai_segments` is kept for WES for symmetry.

The WGS tables came with tx_data PR #11 (2026-11-17). **WES:** the four tables marked "not in tx_data yet" must be
added to tx_data with the same columns as their `wgs_` counterparts before `config.wes.toml` can run (`check`
lists them until then; see `../WGD/TODO.md`, Inputs). **Manual QC** is not a tx_data table: both configs read
the WES manual QC file (`[qc].manual_qc_table`; Piotr, 2026-11-17: every WGS tumour is in the WES cohort).
tx_data's `nemo_mock/` has a minimal mock of it at the mirrored path.

## Sandbox

```bash
cd ../../tx_data && pixi run build sandbox && cd -   # tx_data mock cohort, with tx_data #11
./run.sh config.wgs.sandbox.toml                      # writes output/wgs_sandbox/
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

`tests/make_fixture.py` writes a synthetic store with every table, under either config's table names
(including the proposed WES ones); the tests run for both configs, and the full-run test runs the whole
pipeline on it, PGDD included (skipped without pixi + micromamba). 14 tests (2026-11-24).

## Single-clone tumours

Step 00 sets single-clone tumours aside (they don't go to PGDD). Since 2026-11-17 (Piotr) step 04 reads
them back (`00_pgdd_input/single_clone_samples.tsv`): triage `single_clone`, no events, and the ploidy
fill places GDs on their only clone until every sample matches its `num_gds`, with the ploidy prior,
quality factors and QC as for any tumour. Class `clonal_wgd` or `no_wgd`. Samples that disagree and
can't be reconciled leave the tumour `unresolved`.

## Known gaps / TODO

- WES run blocked on the four WES tx_data tables above.
- TODO (tx_data): give the manual QC mock the real file's full column set.
- `check` stops when `[pgdd].src` is a PGDD without the bootstrap/stability arguments (e.g. the wlippa
  fork at 1af16a8); use amf71/ParallelGDDetect at cb0437b or later.

## Output: quality_pass

`wgd_calls.tsv` keeps every call. `quality_pass` (2026-11-19) is False on a clone carrying a GD with event
quality below `[wgd].quality_cutoff` (the GDs `wgd_calls_rejected_below_<Q>.tsv` drops; mutation-supported GDs
pass unless `reject_mutation_supported`), True otherwise, null on rows without calls or with cutoff -1.
