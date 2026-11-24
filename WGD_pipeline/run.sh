#!/usr/bin/env bash
# WGD_pipeline: PGDD -> WGD caller (WGD steps 00-07) -> ALPACA post-hoc correction, TRACERx WGS or WES.
# Every setting comes from one config file: config.wgs.toml (WGS) or config.wes.toml (WES). Inputs come from
# the tx_data store; [tx_data.tables] picks the cohort's tables.
#
#   ./run.sh <config.toml> [--from STEP] [--to STEP] [--only STEP]
#
# Steps, in order:
#   check    validate the config and list tx_data tables that are missing (stops if any)
#   stage    tx_data tables -> <out>/inputs (layouts the steps read); settings files -> <out>/config
#   wgd00    WGD step 00: per-sample num_gds, PGDD input table (single-clone tumours set aside)
#   pgdd     ParallelGDDetect (detect_par_gd) with [pgdd.settings]
#   wgd01    parse PGDD          wgd02  prune to the ALPACA tree     wgd03  per-sample summary
#   wgd04    place GDs ([wgd.params], QC if [qc].enabled); single-clone tumours by ploidy on their one clone
#   wgd05    per-clone table; GDs below [wgd].quality_cutoff rejected (wgd_calls_rejected_below_<Q>.tsv)
#   wgd07    cohort report (07_report/report.html)
#   correct  ALPACA post-hoc correction on the (filtered) calls
#
# Environment: PIXI (default `pixi`) runs the WGD and correction environments; R comes from
# [pgdd].rscript. Logs: <out>/logs/<step>.log.
set -euo pipefail

usage() { sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }
[[ $# -ge 1 ]] || usage 1
[[ $1 == -h || $1 == --help ]] && usage 0

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd "$HERE/.." && pwd)
CONFIG=$(cd "$(dirname "$1")" && pwd)/$(basename "$1"); shift
STEPS=(check stage wgd00 pgdd wgd01 wgd02 wgd03 wgd04 wgd05 wgd07 correct)
FROM=check; TO=correct
while [[ $# -gt 0 ]]; do
  case $1 in
    --from) FROM=$2; shift 2 ;;
    --to) TO=$2; shift 2 ;;
    --only) FROM=$2; TO=$2; shift 2 ;;
    -h|--help) usage 0 ;;
    *) echo "unknown option: $1" >&2; usage 1 ;;
  esac
done
idx() { local i; for i in "${!STEPS[@]}"; do [[ ${STEPS[$i]} == "$1" ]] && { echo "$i"; return; }; done
        echo "unknown step: $1 (one of: ${STEPS[*]})" >&2; exit 1; }
I_FROM=$(idx "$FROM"); I_TO=$(idx "$TO")

PIXI=${PIXI:-pixi}
WGD="$ROOT/WGD"; APC="$ROOT/alpaca_post_hoc_correction"
wgdpy() { "$PIXI" run --manifest-path "$WGD/pixi.toml" python "$@"; }
apcpy() { "$PIXI" run --manifest-path "$APC/pixi.toml" python "$@"; }
cfg() { wgdpy "$HERE/pipeline/prep.py" get "$CONFIG" "$1"; }

OUT=$(cfg run.out_dir)
COHORT=$(cfg run.cohort); COHORT_LABEL=${COHORT^^}
mkdir -p "$OUT/logs"
IN="$OUT/inputs"; W="$OUT/wgd"; P="$OUT/pgdd"; TX="$IN/tx_data"

run_step() {  # run_step <name> <command...>: log to <out>/logs/<name>.log, stop on failure
  local name=$1; shift
  echo "== $name  $(date -u +%FT%TZ)"
  if ! "$@" > >(tee "$OUT/logs/$name.log") 2>&1; then
    echo "!! step $name failed; see $OUT/logs/$name.log" >&2; exit 1
  fi
}

step_check() { wgdpy "$HERE/pipeline/prep.py" check "$CONFIG"; }
step_stage() { wgdpy "$HERE/pipeline/prep.py" stage "$CONFIG" "$OUT" && wgdpy "$HERE/pipeline/prep.py" write "$CONFIG" "$OUT"; }
step_wgd00() {
  local extra=(); [[ $(cfg wgd.keep_nontree_clusters) == true ]] && extra+=(--keep-nontree-clusters)
  wgdpy "$WGD/scripts/00_prepare_pgdd_input.py" --conipher-dir "$IN/conipher" --msai-dir "$IN/msai" --out "$W" "${extra[@]}"
}
step_pgdd() {
  local rscript; read -r -a rscript <<< "$(cfg pgdd.rscript)"
  "${rscript[@]}" "$HERE/pipeline/run_pgdd.R" "$(cfg pgdd.src)" "$W/00_pgdd_input/pgdd_input_all_tumours.tsv" \
    "$OUT/config/pgdd_settings.json" "$P"
}
step_wgd01() { wgdpy "$WGD/scripts/01_parse_pgdd.py" --pgdd-dir "$P" --out "$W" --expected-settings "$OUT/config/pgdd_settings.json"; }
step_wgd02() { wgdpy "$WGD/scripts/02_prune_tree.py" --tx-data-dir "$TX" --out "$W"; }
step_wgd03() {
  wgdpy "$WGD/scripts/03_generate_summary.py" --pgdd-input "$W/00_pgdd_input/pgdd_input_all_tumours.tsv" \
    --pgdd-regions "$P/pgdd_GDs_per_region.tsv" --out "$W"
}
step_wgd04() {
  local qc=(); [[ $(cfg qc.enabled) == true ]] && qc=(--qc-table "$IN/qc/manual_qc.tsv" --release-dir "$IN/qc/release")
  wgdpy "$WGD/scripts/04_place_events.py" --tx-data-dir "$TX" --out "$W" --alpaca-score "$(cfg wgd.alpaca_score)" \
    --params "$OUT/config/wgd_params.toml" "${qc[@]}"
}
step_wgd05() {
  local q extra=(); q=$(cfg wgd.quality_cutoff)
  awk -v q="$q" 'BEGIN { exit !(q >= 0) }' && extra+=(--reject-below "$q")
  [[ $(cfg wgd.reject_mutation_supported) == true ]] && extra+=(--reject-mutation-supported)
  wgdpy "$WGD/scripts/05_final_table.py" --tx-data-dir "$TX" --out "$W" "${extra[@]}"
}
step_wgd07() {
  wgdpy "$WGD/scripts/07_report.py" --out "$W" --review-quality "$(cfg wgd.review_quality)" --title "WGD calls: TRACERx $COHORT_LABEL"
}
step_correct() {
  local mq extra=(); mq=$(cfg correction.min_quality)
  awk -v q="$mq" 'BEGIN { exit !(q >= 0) }' && extra+=(--min-quality "$mq")
  apcpy "$APC/scripts/correct_alpaca.py" --alpaca "$TX/alpaca.parquet" --clone-proportions "$TX/clone_proportions.parquet" \
    --tree "$TX/alpaca_wgd_ratios.parquet" --wgd-calls "$W/05_table/$(cfg derived.calls_file)" \
    --out "$OUT/correction" "${extra[@]}"
}

echo "WGD_pipeline  cohort=$COHORT  config=$CONFIG  out=$OUT  steps=${STEPS[*]:$I_FROM:$((I_TO - I_FROM + 1))}"
cp "$CONFIG" "$OUT/logs/config.used.toml"
git -C "$ROOT" rev-parse --short HEAD > "$OUT/logs/darwins_analysis_commit.txt" 2>/dev/null || true
for i in $(seq "$I_FROM" "$I_TO"); do run_step "${STEPS[$i]}" "step_${STEPS[$i]}"; done
echo "done: $OUT"
