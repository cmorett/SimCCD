#!/usr/bin/env python3
"""
Lightweight regression checker for paper outputs.

Usage:
  python analysis/regression_check.py --summary paper_outputs/run/tables/validation_summary.csv
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, Tuple

try:
    import pandas as pd
except ImportError:
    pd = None


def load_summary(path: Path) -> Dict[str, float]:
    if pd:
        df = pd.read_csv(path)
        if df.empty:
            return {}
        return df.iloc[0].to_dict()
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            parsed: Dict[str, object] = {}
            for k, v in row.items():
                if v in ("", None):
                    parsed[k] = 0.0
                    continue
                try:
                    parsed[k] = float(v)
                except ValueError:
                    parsed[k] = v
            return parsed
    return {}


def check_range(name: str, value: float, bounds: Tuple[float, float]) -> Tuple[bool, str]:
    lo, hi = bounds
    if value < lo or value > hi:
        return False, f"{name}={value:.3g} outside [{lo}, {hi}]"
    return True, f"{name}={value:.3g}"


def main():
    parser = argparse.ArgumentParser(description="Regression sanity checker for validation_summary.csv.")
    parser.add_argument("--summary", required=True, help="Path to validation_summary.csv")
    parser.add_argument("--min-quality-events", type=int, default=20, help="Minimum non-truncated pixelized events expected")
    parser.add_argument("--max-truncated-fraction", type=float, default=0.1, help="Maximum allowed truncated fraction")
    parser.add_argument("--dedx-range", type=float, nargs=2, default=(0.01, 1e4), help="Allowed dE/dx median range (MeV/cm)")
    parser.add_argument("--sigma-trans-range", type=float, nargs=2, default=(0.05, 300.0), help="Allowed sigma_trans median range (pix)")
    parser.add_argument("--min-hit-fraction", type=float, default=0.0, help="Minimum expected CCD hit fraction (n_events_ccd/n_events)")
    parser.add_argument("--fail-on-warning", action="store_true", help="Exit non-zero if any check fails")
    args = parser.parse_args()

    summary_path = Path(args.summary)
    if not summary_path.exists():
        raise SystemExit(f"Summary file not found: {summary_path}")

    summary = load_summary(summary_path)
    if not summary:
        raise SystemExit(f"No rows found in {summary_path}")

    failures = []
    notes = []

    quality_events = int(summary.get("n_pixel_metrics_quality", 0))
    if quality_events < args.min_quality_events:
        failures.append(f"quality events {quality_events} < {args.min_quality_events}")
    else:
        notes.append(f"quality events = {quality_events}")

    hit_frac = float(summary.get("hit_fraction", 0.0))
    if hit_frac < args.min_hit_fraction:
        failures.append(f"hit fraction {hit_frac:.4f} < {args.min_hit_fraction}")
    else:
        notes.append(f"hit fraction = {hit_frac:.4f}")

    trunc_frac = float(summary.get("truncated_fraction", 0.0))
    ok, msg = check_range("truncated_fraction", trunc_frac, (0.0, args.max_truncated_fraction))
    (notes if ok else failures).append(msg)

    dedx_median = summary.get("dEdx_median", summary.get("dEdx_p50"))
    if dedx_median is not None:
        ok, msg = check_range("dEdx_median", float(dedx_median), tuple(args.dedx_range))
        (notes if ok else failures).append(msg)

    sigma_trans_median = summary.get("sigma_trans_median", summary.get("sigma_trans_p50"))
    if sigma_trans_median is not None:
        ok, msg = check_range("sigma_trans_median", float(sigma_trans_median), tuple(args.sigma_trans_range))
        (notes if ok else failures).append(msg)

    if failures:
        print("Regression check: FAIL")
        for msg in failures:
            print(f"- {msg}")
        for msg in notes:
            print(f"* {msg}")
        if args.fail_on_warning:
            sys.exit(1)
    else:
        print("Regression check: PASS")
        for msg in notes:
            print(f"- {msg}")


if __name__ == "__main__":
    main()
