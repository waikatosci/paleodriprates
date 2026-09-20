"""
make_canonical_te_input.py — build the canonical Dr Paleo TE input from the
raw HS4 table by applying the documented exclusion list.

Decision record (2026-07-20): the canonical input is the FULL raw table
(dr_app/HS4_example_inputs/HS4_TE.csv, 589 rows, value-identical to
Drip_rate.xlsx sheet 2.Trace_Elems) minus the rows listed in
excluded_points.csv. As of the decision date that is a single exclusion —
the 0.06 cm surface-cap point — giving 588 rows spanning 0.18–255.5 cm.

The legacy 13-row exclusion behind the published 576-point external grid is
RETIRED: twelve of those rows are not reproducible by the pipeline's own
Hampel/Pareto detector and had no recoverable rationale, so they are
re-included. Spikes among them are handled by the run's detector, which
median-replaces (never drops): on the full series it flags Ni @ 111.10 cm
and Co @ 8.48, 36.03 cm at production settings (window 11, ZERO_TOL 1e-3,
PDF_TOL 3). Note for the record: all twelve re-included depths sit above
their local Ni baseline (1.24-2.59x) — see the session log of 2026-07-20;
if a detrital-contamination screen is ever re-instated, add rows to
excluded_points.csv and re-run this script rather than editing data files.

The MS statement "n = 585" was confirmed an error (not a third cleaning
rule); restate n from this pipeline at assembly: 588 on the depth grid,
586 within the dated span (age model 0.0–253.0 cm excludes the rows at
253.5 and 255.5 for age-mode products).

USAGE
    python3 make_canonical_te_input.py \
        [--raw ../dr_app/HS4_example_inputs/HS4_TE.csv] \
        [--exclusions excluded_points.csv] \
        [--out ../dr_app/HS4_example_inputs/HS4_TE_canonical.csv]
"""
import argparse
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default=os.path.join(
        HERE, "..", "dr_app", "HS4_example_inputs", "HS4_TE.csv"))
    ap.add_argument("--exclusions", default=os.path.join(
        HERE, "excluded_points.csv"))
    ap.add_argument("--out", default=os.path.join(
        HERE, "..", "dr_app", "HS4_example_inputs", "HS4_TE_canonical.csv"))
    args = ap.parse_args()

    raw = pd.read_csv(args.raw)
    exc = pd.read_csv(args.exclusions)

    depth_col, ni_col = raw.columns[0], raw.columns[1]
    depth = np.round(raw[depth_col].astype(float), 4)
    ni = raw[ni_col].astype(float)
    # An exclusion matches on depth; where the list gives ni_ppm as well, the
    # Ni value must also match (within 0.01 ppm). This disambiguates the three
    # depths at which the raw table carries one row from each laboratory
    # (8.09, 156.38, 157.48 cm) so that only the second-laboratory row is dropped.
    keep = np.ones(len(raw), dtype=bool)
    for _, e in exc.iterrows():
        m = depth == round(float(e["depth_cm"]), 4)
        if "ni_ppm" in exc.columns and pd.notna(e.get("ni_ppm")):
            m &= np.isclose(ni, float(e["ni_ppm"]), atol=0.01)
        keep &= ~m
    out = raw[keep].reset_index(drop=True)

    n_dropped = int((~keep).sum())
    if n_dropped != len(exc):
        raise SystemExit(
            f"exclusion mismatch: {len(exc)} rows listed, {n_dropped} matched "
            f"in {args.raw} — check depths.")

    out.to_csv(args.out, index=False)
    print(f"raw rows: {len(raw)}  excluded: {n_dropped}  canonical: {len(out)}")
    print(f"depth span: {out[depth_col].min()}–{out[depth_col].max()} cm")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
