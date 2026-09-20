#!/usr/bin/env python3
"""
export_multielement_te.py — export the full multi-element stalagmite TE table
(V, Cr, Co, Ni, Cu, Zn in calcite, ppm) from the canonical workbook sheet
HS4_SourceData.xlsx[01_master_TE] to a flat CSV for the companion screens.

Applies the same documented exclusions as the drip-rate canonical input:
the 0.06 cm surface-cap point (calibration/excluded_points.csv) and the 20
second-laboratory samples flagged on 01_master_TE column "excluded"
(2026-09-21), so the screened series and the inverted series describe the
same sample set. Values <= 0 (below detection; occurs for Zn) are exported as
blank. Rows are sorted by depth; duplicate depths (parallel drill/mill
transects at 8.09, 156.38, 157.48 cm) are retained, as in the canonical
input.

Output: ../manuscript_figures/external/HS4_TE_multielement.csv
Columns: sample, depth_cm, age_yBP, Ni, Co, Cu, Cr, V, Zn   (calcite ppm)
"""
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
WB = os.path.join(ROOT, "manuscript_figures", "HS4_SourceData.xlsx")
EXC = os.path.join(HERE, "excluded_points.csv")
OUT = os.path.join(ROOT, "manuscript_figures", "external",
                   "HS4_TE_multielement.csv")

COLMAP = {"Ni_calcite_ppm": "Ni", "Co_calcite_ppm": "Co",
          "Cu_calcite_ppm": "Cu", "Cr_calcite_ppm": "Cr",
          "V_calcite_ppm": "V", "Zn_calcite_ppm": "Zn"}


def main():
    te = pd.read_excel(WB, sheet_name="01_master_TE")
    te = te.dropna(subset=["depth_cm"])
    exc = pd.read_csv(EXC)
    # depth-only exclusions (the 0.06 cm surface-cap point); the second-
    # laboratory rows are flagged on the sheet itself (column "excluded")
    # because three of them share a depth with a primary-laboratory row.
    depth_only = exc[exc.get("ni_ppm", pd.Series(dtype=float)).isna()] if "ni_ppm" in exc.columns else exc
    drop = set(np.round(depth_only["depth_cm"].astype(float), 4))
    te = te[~np.round(te.depth_cm.astype(float), 4).isin(drop)]
    if "excluded" in te.columns:
        n_flag = te["excluded"].notna().sum()
        te = te[te["excluded"].isna()]
        print(f"dropped {n_flag} rows flagged on 01_master_TE[excluded]")

    out = te[["sample", "depth_cm", "age_yBP_new"] + list(COLMAP)].copy()
    out = out.rename(columns={"age_yBP_new": "age_yBP", **COLMAP})
    for e in COLMAP.values():
        out.loc[out[e] <= 0, e] = np.nan     # below detection
    out = out.sort_values("depth_cm").reset_index(drop=True)
    out.to_csv(OUT, index=False, float_format="%.6g")
    print(f"wrote {OUT}: {len(out)} rows, depth "
          f"{out.depth_cm.min():.2f}-{out.depth_cm.max():.2f} cm, "
          f"{out.age_yBP.notna().sum()} dated")


if __name__ == "__main__":
    main()
