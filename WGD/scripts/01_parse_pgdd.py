"""Step 01: read PGDD outputs and write long-format parquet tables.

Usage:
  pixi run python scripts/01_parse_pgdd.py --nemo-root /nemo --out output/<name>
  pixi run python scripts/01_parse_pgdd.py --pgdd-dir /path/to/pgdd --out output/<name>

`--nemo-root` defaults to $TX_NEMO_ROOT. Locally, point it at tx_data's
`nemo_mock/` (needs tx_data PR #5, branch darwin/pgdd-mock).
`--expected-settings <json>`: PGDD settings to check the run against (WGD_WGS_pipeline passes its
[pgdd.settings]); default `wgd.pgdd.EXPECTED_SETTINGS`. A `run.json` in the PGDD directory (written by
the pipeline's PGDD wrapper) supplies values the outputs don't carry.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from wgd import pgdd  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group()
    src.add_argument("--nemo-root", type=Path, default=os.environ.get("TX_NEMO_ROOT"))
    src.add_argument("--pgdd-dir", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--expected-settings", type=Path, help="JSON of the PGDD settings the run should have used")
    args = ap.parse_args()

    if args.pgdd_dir is None and args.nemo_root is None:
        ap.error("give --pgdd-dir or --nemo-root (or set TX_NEMO_ROOT)")
    directory = args.pgdd_dir or pgdd.pgdd_dir(args.nemo_root)

    tables = pgdd.read_pgdd(directory)
    out = args.out / "01_pgdd"
    out.mkdir(parents=True, exist_ok=True)

    tables["tumours"].write_parquet(out / "pgdd_tumours.parquet")
    tables["regions"].write_parquet(out / "pgdd_regions.parquet")
    tables["events"].write_parquet(out / "pgdd_events.parquet")
    tables["mut_counts"].write_parquet(out / "pgdd_mut_counts.parquet")
    tables["bootstrap"].write_parquet(out / "pgdd_bootstrap.parquet")
    pgdd.event_clusters(tables["events"]).write_parquet(out / "pgdd_event_clusters.parquet")
    pgdd.event_regions(tables["events"]).write_parquet(out / "pgdd_event_regions.parquet")
    pgdd.region_clusters(tables["regions"]).write_parquet(out / "pgdd_region_clusters.parquet")

    summary = pgdd.summarise(tables)
    expected = json.loads(args.expected_settings.read_text()) if args.expected_settings else None
    rec_path = Path(directory) / "run.json"
    recorded = json.loads(rec_path.read_text()).get("settings") if rec_path.exists() else None
    settings_diff = pgdd.check_settings(tables, expected, recorded)
    run = {
        "pgdd_dir": str(directory),
        "settings": tables["settings"].to_dicts()[0],
        "settings_differing_from_expected": settings_diff,
        "summary": summary,
    }
    (out / "run.json").write_text(json.dumps(run, indent=2, default=str) + "\n")
    for name, d in settings_diff.items():
        print(f"WARNING: PGDD {name} = {d['found']}, expected {d['expected']}")
    for k, v in summary.items():
        print(f"{k:>26}: {v}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
