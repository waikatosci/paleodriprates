#!/usr/bin/env python3
"""
crosslab_sensitivity.py -- Supplementary Figure 17 and the numbers behind
Supplementary Methods 14.3: what the excluded second-laboratory samples would
imply for the 5.2 ka event if they were retained after a cross-laboratory
correction.

Twenty samples in the working trace-element table came from a supplementary
analytical run at a second laboratory (12 with prefix HS4-A- across the 5.2 ka
event, 155.2-157.8 cm; 8 with prefix HS4-C- at the core top, 7.0-8.5 cm). The
data owner reports that the second-laboratory calibration does not reproduce a
certified standard (Ni 5.20 ppm returned 2.38 ppm) and asked for the samples
to be removed; the primary reconstruction therefore uses the 568-point
single-laboratory input (calibration/excluded_points.csv).

At three depths both laboratories analysed the same horizon (8.09, 156.38,
157.48 cm). Across those pairs the second-laboratory values are a linear
function of the primary-laboratory values (r >= 0.998 for Ni and Co),
consistent with a near-constant additive offset (~+3.9 ppm Ni, ~+0.9 ppm Co).
Rescaling the second-laboratory rows with that fit and re-running the
inversion gives an upper-bound sensitivity on the event amplitude. The fit is
not independently validated (the certified-standard test points the other way,
an under-read), so it is reported as a bound, not adopted.

Inputs   dr_app/HS4_example_inputs/HS4_TE.csv          raw table (both laboratories)
         calibration/excluded_points.csv               the 20 excluded rows
         runs written by drive_run.py: hr_clean (568-pt) and hr_corr (585-pt, rescaled)
Outputs  ../manuscript_figures/output/FigS_crosslab_sensitivity.{png,pdf}
         ../manuscript_figures/output/TableS_crosslab_pairs.csv
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(ROOT, "manuscript_figures"))
from ngeo_style import (apply_style, style_ax, panel_label, label_bubble,  # noqa: E402
                        COL_NI, COL_CO, COL_BG_900, COL_BG_600, COL_BG_300,
                        COL_BG_50, COL_DORANGE, COL_TEAL_600, DOUBLE_COL)
import matplotlib.pyplot as plt  # noqa: E402

RAW = os.path.join(ROOT, "dr_app", "HS4_example_inputs", "HS4_TE.csv")
EXC = os.path.join(ROOT, "calibration", "excluded_points.csv")
RUNS = os.environ.get("DRPALEO_RUNS", "/home/claude/work/runs")
OUT = os.path.join(ROOT, "manuscript_figures", "output")
PAIR_DEPTHS = [8.09, 156.38, 157.48]
MODERN = 16.66


def main():
    os.makedirs(OUT, exist_ok=True)
    raw = pd.read_csv(RAW)
    d, ni, co = raw.columns[:3]
    exc = pd.read_csv(EXC)
    exc = exc[exc["ni_ppm"].notna()]
    second = np.zeros(len(raw), bool)
    for _, e in exc.iterrows():
        second |= (np.round(raw[d], 4) == round(float(e.depth_cm), 4)) & np.isclose(raw[ni], e.ni_ppm, atol=0.01)
    S, P = raw[second], raw[~second]

    # ── paired depths and the cross-laboratory fit ───────────────────────
    pairs = []
    for dep in PAIR_DEPTHS:
        a = P[np.isclose(P[d], dep, atol=0.006)].iloc[0]
        s = S[np.isclose(S[d], dep, atol=0.006)].iloc[0]
        pairs.append(dict(depth_cm=dep, Ni_primary=a[ni], Co_primary=a[co],
                          Ni_second=s[ni], Co_second=s[co]))
    pairs = pd.DataFrame(pairs)
    fit = {}
    for el in ("Ni", "Co"):
        x, y = pairs[f"{el}_second"], pairs[f"{el}_primary"]
        a, b = np.polyfit(x, y, 1)
        r = np.corrcoef(x, y)[0, 1]
        fit[el] = (a, b, r)
        print(f"{el}: primary = {a:.4f} x second {b:+.4f}   r = {r:.4f}   "
              f"(mean second/primary ratio {np.mean(x / y):.2f}; mean offset {np.mean(x - y):+.2f} ppm)")
    pairs.to_csv(os.path.join(OUT, "TableS_crosslab_pairs.csv"), index=False)

    # ── inversions: cleaned (568) vs corrected-inclusive (585) ───────────
    clean = pd.read_csv(os.path.join(RUNS, "hr_clean", "drip_rate_summary.csv"))
    corr = pd.read_csv(os.path.join(RUNS, "hr_corr", "drip_rate_summary.csv"))
    ref = lambda df: df[((df.depth >= 140) & (df.depth < 155)) | ((df.depth > 158) & (df.depth <= 175))].pc50.median()
    for lab, df in (("cleaned 568", clean), ("corrected 585", corr)):
        w = df[(df.depth >= 150) & (df.depth <= 165)]
        m = w.loc[w.pc50.idxmin()]
        print(f"{lab}: minimum {m.pc50:.2f} drips/min at {m.depth:.2f} cm "
              f"({100 * (1 - m.pc50 / ref(df)):.0f} % below the surrounding mid-Holocene median {ref(df):.1f})")

    # ── figure ───────────────────────────────────────────────────────────
    apply_style()
    fig, axs = plt.subplots(1, 3, figsize=(DOUBLE_COL, DOUBLE_COL * 0.34))

    # (a) paired-depth relation
    ax = axs[0]
    for el, col, mk in (("Ni", COL_NI, "o"), ("Co", COL_CO, "s")):
        x, y = pairs[f"{el}_second"], pairs[f"{el}_primary"]
        a, b, r = fit[el]
        xx = np.linspace(0, 9, 50)
        ax.plot(xx, a * xx + b, color=col, lw=0.8)
        ax.plot(x, y, mk, ms=3.5, color=col, label=f"{el} (r = {r:.3f})")
    ax.plot([0, 9], [0, 9], color=COL_BG_300, lw=0.5, ls=":")
    ax.set_xlim(0, 9); ax.set_ylim(0, 6)
    ax.set_xlabel("Second laboratory (ppm)")
    ax.set_ylabel("Primary laboratory (ppm)")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2, fontsize=5.2, handlelength=1.4)
    label_bubble(ax, "1:1", xy=(5.6, 5.6), xytext=(7.4, 4.3), va="center",
                 color=COL_BG_600, bubble_ec=COL_BG_600, arrow_color=COL_BG_600)
    panel_label(ax, "a")

    # (b) concentrations across the 5.2 ka window
    ax = axs[1]
    win = (raw[d] >= 153.5) & (raw[d] <= 159.5)
    Pw, Sw = P[win[~second]], S[win[second]]
    for el, col in (("Ni", COL_NI), ("Co", COL_CO)):
        c = ni if el == "Ni" else co
        a, b, _ = fit[el]
        ax.plot(Pw[d], Pw[c], "o-", ms=2.6, lw=0.6, color=col, label=f"{el}, primary")
        ax.plot(Sw[d], Sw[c], "x", ms=3.2, mew=0.7, color=col, alpha=0.6, label=f"{el}, second (as reported)")
        ax.plot(Sw[d], a * Sw[c] + b, "D", ms=2.4, mfc="none", mew=0.7, color=col, label=f"{el}, second (rescaled)")
    ax.set_xlim(153.5, 159.5); ax.set_yscale("log"); ax.set_ylim(0.2, 15)
    ax.set_xlabel("Depth (cm)")
    ax.set_ylabel("Concentration in calcite (ppm)")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3, fontsize=4.8, columnspacing=0.8, handlelength=1.4)
    panel_label(ax, "b")

    # (c) inverted drip rate, cleaned vs corrected-inclusive
    ax = axs[2]
    for lab, df, col, mk in (("568-point record (used)", clean, COL_BG_900, "o"),
                             ("second laboratory retained, rescaled", corr, COL_DORANGE, "D")):
        w = df[(df.depth >= 153.5) & (df.depth <= 159.5)].sort_values("depth")
        ax.fill_between(w.depth, w.pc25, w.pc75, color=col, alpha=0.15, lw=0)
        ax.plot(w.depth, w.pc50, marker=mk, ms=2.4, lw=0.7, color=col, mfc=("none" if mk == "D" else col), label=lab)
    ax.axhline(MODERN, color=COL_BG_300, lw=0.5, ls="--")
    label_bubble(ax, "modern baseflow", xy=(158.6, MODERN), xytext=(157.9, 40), va="center",
                 color=COL_BG_600, bubble_ec=COL_BG_600, arrow_color=COL_BG_600)
    ax.set_xlim(153.5, 159.5); ax.set_yscale("log"); ax.set_ylim(1, 60)
    ax.set_xlabel("Depth (cm)")
    ax.set_ylabel("Inferred drip rate (drips min$^{-1}$)")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=1, fontsize=5.2, handlelength=1.4)
    panel_label(ax, "c")

    for ax in axs:
        style_ax(ax)
    fig.tight_layout(w_pad=0.6)
    fig.subplots_adjust(bottom=0.24, left=0.055, right=0.985, wspace=0.34)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"FigS_crosslab_sensitivity.{ext}"))
    print(f"wrote {os.path.join(OUT, 'FigS_crosslab_sensitivity.png')} (+pdf)")


if __name__ == "__main__":
    main()
