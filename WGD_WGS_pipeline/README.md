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

| table | used by | in tx_data |
|---|---|---|
| `wgs_consensus_segments` | step 00 (num_gds) | ✅ |
| `wgs_clone_proportions` | steps 02, 04, correction | ✅ |
| `wgs_conipher_mutations` | step 00 (SCoutput.CLEAN + treeTable, per mutation × sample) | ❌ to add |
| `wgs_conipher_tree_edges` | step 00 (default tree) | ❌ to add |
| `wgs_alpaca` | step 04, correction (clone CN per segment; WES `alpaca` columns) | ❌ to add |
| `wgs_alpaca_wgd_ratios` | steps 02, 04, 05, correction (tree + ratios; WES `alpaca_wgd_ratios` columns) | ❌ to add |
| `wgs_manual_qc`, `wgs_conipher_cluster_info` | step 04 QC, only with `[qc].enabled` | ❌ to add |

Proposed schemas for the missing tables: `pipeline/wgs.py` (top).

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
  `[pgdd.settings]`: changing a setting here gives a warning there.
