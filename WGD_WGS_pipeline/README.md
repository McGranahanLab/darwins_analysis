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

On the mock (2026-11-17): 11 tumours, 8 single-clone (set aside, see Known gaps); LTX0004 gets a clonal
GD; LTX0103 a borderline clonal GD (genome fraction 0.53) that the 0.3 cut-off rejects;
LTX0104-Tumour2 is `qc_failed` (truncal CCF 86 in the mock). Step 01 warns that the PGDD settings
are missing because PGDD's `mut_counts` is empty when no cluster is tested; the settings used are
in `pgdd/run.json`.

## Requirements

`pixi` (runs the WGD and correction environments), R with `data.table` + `jsonlite` (`[pgdd].rscript`,
e.g. `micromamba run -n pgdd Rscript`), and an amf71/ParallelGDDetect checkout at `cb0437b` (`[pgdd].src`).

## Tests

```bash
pixi run --manifest-path ../WGD/pixi.toml pytest tests/ -q
```

`tests/make_fixture.py` writes a synthetic WGS store with every table (including the proposed ones);
the last test runs the whole pipeline on it, PGDD included (skipped without pixi + micromamba).

## Known gaps

- **Single-clone tumours** never reach `wgd_calls.tsv`: step 00 sets them aside
  (`00_pgdd_input/single_clone_events.tsv`) and steps 04–05 don't read them back, so the correction
  skips them. This affects WES too, and matters more for WGS.
- Step 01 checks the PGDD run against the settings hard-coded in `WGD/wgd/pgdd.py`, not against
  `[pgdd.settings]`: changing a setting here gives a warning there (and every setting warns when
  PGDD's `mut_counts` is empty).
- TODO (tx_data): give the manual QC mock the real file's full column set.
