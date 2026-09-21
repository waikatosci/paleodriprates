#!/usr/bin/env python3
"""
make_fig_events_differ.py — the "8.2 ka and 5.2 ka differ in kind" figure
(Supplementary Figure 10; Supplementary Methods 13). Realises the four-panel
figure the SI caption specifies.

NCOMMS-26-041445-T. Answers Reviewer 2 major 2 (8.2 ka reconciliation) and
Reviewer 3 major 2 (exclude detrital/PCP at 5.2 ka) by CONTRASTING the two
intervals. At 5.2 ka the kinetic Co-Ni coupling is pristine (record maximum)
and the isotopes are silent; at 8.2 ka the chronology breaks, the Co-Ni
coupling collapses, and the isotopes/PCP indicators excurse.

Panels (all reproducible from the released Source Data + companion CSV):
  (a) U-Th age-depth relation across the base: four stratigraphic inversions
      confined to 228-236.8 cm (10 of the 30 U-Th dates), 2sigma median there 200 yr vs 60 yr elsewhere.
  (b) Co vs Ni within the 8.2 ka interval (flanking samples grey): coupling
      collapsed (Spearman rho on logs ~ 0.05).
  (c) Co vs Ni within the 5.2 ka interval: coupling pristine, the strongest
      in the record (rho ~ 0.96).
  (d) Co-Ni coupling in 0.4 ka sliding windows across the Holocene: 5.2 ka is
      the maximum, 8.2 ka a minimum, against a record median ~ 0.46; coupling
      is uncorrelated with Co concentration (not an amplitude effect).

Inputs   ../manuscript_figures/HS4_SourceData.xlsx  (02_chronology, 03_isotopes)
         ../manuscript_figures/external/HS4_TE_multielement.csv
Outputs  ../manuscript_figures/output/FigS_events_differ.{png,pdf}
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(ROOT, "manuscript_figures"))
from ngeo_style import (panel_label, COL_DORANGE, COL_BROWN_200, COL_BG_200, COL_BG_600, apply_style, style_ax, COL_NI, COL_CO, COL_BG_900,   # noqa: E402
                        COL_BG_300, COL_BG_500, COL_TEAL_600, COL_BROWN_500,
                        COL_RED_900, DOUBLE_COL)

import matplotlib.pyplot as plt                                              # noqa: E402

WB = os.path.join(ROOT, "manuscript_figures", "HS4_SourceData.xlsx")
TE_CSV = os.path.join(ROOT, "manuscript_figures", "external", "HS4_TE_multielement.csv")
OUTDIR = os.path.join(ROOT, "manuscript_figures", "output")

WIN_82 = (230.8, 236.8)     # cm — matches SM13.2
WIN_52 = (155.0, 157.9)     # cm — the 5.2 ka excursion
SLIDE_HALF = 200            # yr (0.4 ka windows)


def load():
    ch = pd.read_excel(WB, sheet_name="02_chronology", header=6)
    tp = ch[ch.is_UTh_tiepoint.astype(str).str.lower().eq("yes")].copy()
    te = pd.read_csv(TE_CSV)
    dated = te.dropna(subset=["age_yBP"]).sort_values("age_yBP")
    return tp, te, dated


def logrho(d):
    return stats.spearmanr(np.log(d.Ni), np.log(d.Co))[0]


def sliding(dated):
    centers = np.arange(300, 8700, 50)
    rows = []
    for c in centers:
        w = dated[(dated.age_yBP >= c - SLIDE_HALF) & (dated.age_yBP <= c + SLIDE_HALF)]
        if len(w) >= 6:
            rows.append((c, logrho(w), w.Co.median()))
    return pd.DataFrame(rows, columns=["age", "rho", "co"])


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    tp, te, dated = load()

    base = tp[(tp.depth_cm >= 226) & (tp.depth_cm <= 238)].sort_values("depth_cm")
    inv = base[base.age_yBP.diff() < 0]
    e82 = te[(te.depth_cm >= WIN_82[0]) & (te.depth_cm <= WIN_82[1])]
    e52 = te[(te.depth_cm >= WIN_52[0]) & (te.depth_cm <= WIN_52[1])]
    sl = sliding(dated)
    r_co, p_co = stats.spearmanr(sl.rho, sl.co)

    print(f"(a) inversions {inv.depth_cm.round(1).tolist()}; 2sigma there "
          f"{base.err2s_yr.median():.0f} vs "
          f"{tp[~tp.index.isin(base.index)].err2s_yr.median():.0f} yr")
    print(f"(b) 8.2 ka coupling rho = {logrho(e82):.2f} (n={len(e82)})")
    print(f"(c) 5.2 ka coupling rho = {logrho(e52):.2f} (n={len(e52)}) "
          f"[record max = {sl.rho.max():.2f}]")
    print(f"(d) sliding median rho = {sl.rho.median():.2f}; 8.2 ka min "
          f"{sl[(sl.age >= 8000) & (sl.age <= 8300)].rho.min():.2f}; "
          f"coupling vs Co-conc rho = {r_co:.2f} (p = {p_co:.2f}) — not amplitude")

    apply_style()
    fig, axs = plt.subplots(2, 2, figsize=(DOUBLE_COL, 4.4))

    # (a) age-depth inversions
    ax = axs[0, 0]
    ax.axvspan(228, 236.8, color=COL_BG_300, alpha=0.25, lw=0)
    a = tp.sort_values("depth_cm")
    ax.errorbar(a.depth_cm, a.age_yBP / 1000, yerr=a.err2s_yr / 1000, fmt="o",
                ms=2.4, lw=0.5, color=COL_BG_900, ecolor=COL_BG_500, capsize=1)
    ax.scatter(inv.depth_cm, inv.age_yBP / 1000, s=32, facecolor="none",
               edgecolor=COL_DORANGE, lw=0.9, zorder=6)
    ax.set_xlim(200, 256)
    ax.set_xlabel("Depth (cm)")
    ax.set_ylabel("Age (ka BP)")
    ax.set_ylim(6.9, 9.8)
    med_in = base.err2s_yr.median(); med_out = tp[~tp.index.isin(base.index)].err2s_yr.median()
    ax.annotate(f"{len(inv)} inversions\n2\u03c3 {med_in:.0f} vs {med_out:.0f} yr", xy=(232, 7.55), fontsize=5, ha="center", va="top", color=COL_BG_600)
    panel_label(ax, "a")

    # (b) Co vs Ni at 8.2 ka
    ax = axs[0, 1]
    fl = dated[(dated.age_yBP >= 7600) & (dated.age_yBP <= 8700)]
    fl = fl[~fl.index.isin(e82.index)]
    ax.scatter(fl.Ni, fl.Co, s=6, color=COL_BG_300, alpha=0.6, lw=0, label="flanking")
    ax.scatter(e82.Ni, e82.Co, s=14, color=COL_BROWN_500, lw=0.3,
               edgecolor=COL_BG_900, label="8.2 ka", zorder=5)
    ax.text(0.10, 0.90, f"\u03c1 = {logrho(e82):+.2f}", transform=ax.transAxes,
            fontsize=6, color=COL_BROWN_500)
    ax.set_xlabel("Ni (ppm)")
    ax.set_ylabel("Co (ppm)")
    ax.legend(frameon=False, fontsize=5, loc="lower right")
    panel_label(ax, "b")

    # (c) Co vs Ni at 5.2 ka
    ax = axs[1, 0]
    fl2 = dated[(dated.age_yBP >= 4700) & (dated.age_yBP <= 5600)]
    fl2 = fl2[~fl2.index.isin(e52.index)]
    ax.scatter(fl2.Ni, fl2.Co, s=6, color=COL_BG_300, alpha=0.6, lw=0, label="flanking")
    ax.scatter(e52.Ni, e52.Co, s=14, color=COL_TEAL_600, lw=0.3,
               edgecolor=COL_BG_900, label="5.2 ka", zorder=5)
    b1, b0 = np.polyfit(e52.Ni, e52.Co, 1)
    xs = np.linspace(e52.Ni.min(), e52.Ni.max(), 10)
    ax.plot(xs, b0 + b1 * xs, color=COL_BG_900, lw=0.6, ls="--")
    ax.text(0.10, 0.90, f"\u03c1 = {logrho(e52):+.2f}", transform=ax.transAxes,
            fontsize=6, color=COL_TEAL_600)
    ax.set_xlabel("Ni (ppm)")
    ax.set_ylabel("Co (ppm)")
    ax.legend(frameon=False, fontsize=5, loc="lower right")
    panel_label(ax, "c")

    # (d) sliding coupling
    ax = axs[1, 1]
    ax.plot(sl.age / 1000, sl.rho, color=COL_BG_900, lw=0.8)
    ax.axhline(sl.rho.median(), color=COL_BG_500, lw=0.5, ls=":")
    ax.annotate(f"median {sl.rho.median():.2f}", xy=(0.3, sl.rho.median() + 0.02),
                fontsize=5, color=COL_BG_600, ha="right")
    for lo, hi, c in ((5.0, 5.2, COL_BROWN_200), (8.0, 8.3, COL_BG_200)):
        ax.axvspan(lo, hi, color=c, alpha=0.6, lw=0, zorder=0.5)
    imax = sl.rho.idxmax()
    ax.scatter([sl.loc[imax, "age"] / 1000], [sl.rho.max()], s=16,
               color=COL_TEAL_600, zorder=5)
    ax.set_xlim(9, 0)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Age (ka BP)")
    ax.set_ylabel("Co–Ni coupling ρ (0.4 ka window)")
    panel_label(ax, "d")

    for ax in axs.flat:
        style_ax(ax)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUTDIR, f"FigS_events_differ.{ext}"))
    print(f"\nwrote {os.path.join(OUTDIR, 'FigS_events_differ.png')} (+pdf)")


if __name__ == "__main__":
    main()
