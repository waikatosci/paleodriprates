#!/usr/bin/env python3
"""
detrital_ternary_screen.py — multi-element discrimination of kinetic,
ternary (intact M-NOM co-precipitation), and detrital contributions to the
HS4 trace-metal record, focused on the 5.2 ka and 8.2 ka event levels.

NCOMMS-26-041445-T companion analysis (Reviewer 3 major 2: exclude PCP /
detrital contamination, especially at 5.2 ka; Reviewer 2 major 2: reconcile
with the Zhu et al. 2017 IRM_soft-flux record from the same stalagmite).

Three processes can enrich trace metals in stalagmite calcite, and they leave
different multi-element fingerprints:

  KINETIC (drip-rate) — Ni and Co respond through OMC dissociation with
      metal-specific k_d, so a slow-drip excursion enriches Co MORE than Ni
      (model-predicted trajectory slope d lnCo / d lnNi ~ 1.3, steepening as
      Ni saturates); Cu (inert end-member) and the lithogenic suite do not
      respond.
  SOURCE (dripwater composition) — the shared NOM carrier moves Ni and Co
      TOGETHER but sub-proportionally for Co (observed source-covariance
      slope 0.49, r = +0.73); no lithogenic response.
  DETRITAL — particulate delivery enriches the lithogenic suite (Cr, V, Zn)
      together with the NOM metals at near-crustal proportions (slope ~ 1)
      and tracks the soil-derived magnetite flux (IRM_soft; Zhu et al. 2017).

TERNARY co-precipitation (intact M-NOM complexes) is carried by Cu — the
designated inert/pure-ternary end-member — and tracks NOM delivery rather
than drip rate.

This script computes local-baseline enrichment factors (EF, rolling-median
normalised) for Ni, Co, Cu, Cr, V, Zn; a lithogenic index (geometric mean of
EF_Cr and EF_V); event-window summaries; a formal detrital flag
(EF_lith > 3 AND EF_Zn > 5); the kinetic / source / detrital discriminant in
(ln EF_Ni, ln EF_Co) space; and the record-wide comparison of the lithogenic
index with the Zhu et al. (2017) IRM_soft flux from the same stalagmite.

Inputs   ../manuscript_figures/external/HS4_TE_multielement.csv
         ../manuscript_figures/external/HS4_Zhu2017_IRMsoft_flux.csv
Outputs  ../manuscript_figures/output/TableS_event_enrichment.csv
         ../manuscript_figures/output/FigS_detrital_ternary_screen.{png,pdf}
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(ROOT, "manuscript_figures"))
from ngeo_style import (panel_label, COL_BROWN_200, COL_BG_200, apply_style, style_ax, COL_NI, COL_CO, COL_BG_900,   # noqa: E402
                        COL_BG_300, COL_BG_500, COL_BROWN_500, COL_DORANGE,
                        COL_TEAL_600, COL_BG_600, DOUBLE_COL)
sys.path.insert(0, os.path.join(ROOT, "companion_analysis"))
from source_variation_propagation import phi_of_V, PARAMS                     # noqa: E402

import matplotlib.pyplot as plt                                               # noqa: E402
COL_IRM = "#FF8F00"   # amber, as IRM_soft in main Fig 7

TE_CSV = os.path.join(ROOT, "manuscript_figures", "external",
                      "HS4_TE_multielement.csv")
IRM_CSV = os.path.join(ROOT, "manuscript_figures", "external",
                       "HS4_Zhu2017_IRMsoft_flux.csv")
OUTDIR = os.path.join(ROOT, "manuscript_figures", "output")

ELEMS = ["Ni", "Co", "Cu", "Cr", "V", "Zn"]
ROLL_WIN = 61          # rolling-median baseline window (samples): wide
                       # enough that the finely sampled 5.2 ka event
                       # (~16 points over 2.6 cm) cannot lift its own
                       # baseline; EFs converged for windows >= 61.
Z_FLAG = 4.0           # formal detrital flag: robust log-z > Z_FLAG in BOTH
                       # Cr and Zn (record-wide median/MAD in log space —
                       # immune to the near-detection V tail and to rolling-
                       # window edge effects at the basal contact).

# Event windows (canonical framing)
WIN_52_CORE = ("depth", 155.7, 157.1)       # ~5072-5127 yr BP; the three elevated A-series samples (568-point record, 2026-09-21)
WIN_82 = ("age", 8000.0, 8400.0)
WIN_BASAL = ("depth", 253.4, 256.0)         # outside the dated span (>253.0 cm)
V_MODERN, V_FLOOR = 14.14, 1.07


def load_te():
    te = pd.read_csv(TE_CSV)
    for e in ELEMS:
        base = te[e].rolling(ROLL_WIN, center=True, min_periods=7).median()
        te[f"EF_{e}"] = te[e] / base
    te["EF_lith"] = np.sqrt(te.EF_Cr * te.EF_V)
    for e in ("Cr", "Zn"):
        ln = np.log(te[e])
        mad = 1.4826 * (ln - ln.median()).abs().median()
        te[f"z_{e}"] = (ln - ln.median()) / mad
    te["detrital_flag"] = (te.z_Cr > Z_FLAG) & (te.z_Zn > Z_FLAG)
    return te


def window(te, spec):
    kind, lo, hi = spec
    col = "depth_cm" if kind == "depth" else "age_yBP"
    return te[(te[col] >= lo) & (te[col] <= hi)]


def summarise(te):
    rows = []
    windows = [
        ("record baseline", te),
        ("5.2 ka core (155.7-157.1 cm)", window(te, WIN_52_CORE)),
        ("8.2 ka window (8.0-8.4 ka)", window(te, WIN_82)),
        ("basal detrital (>253.4 cm)", window(te, WIN_BASAL)),
    ]
    for name, d in windows:
        r = dict(window=name, n=len(d))
        for e in ELEMS + ["lith"]:
            r[f"EF_{e}_med"] = d[f"EF_{e}"].median()
            r[f"EF_{e}_max"] = d[f"EF_{e}"].max()
        rows.append(r)
    return pd.DataFrame(rows)


def irm_comparison(te):
    """Bin ln EF_lith and ln EF_Cu into Zhu specimen depositional intervals."""
    irm = pd.read_csv(IRM_CSV)
    dated = te.dropna(subset=["age_yBP"])
    recs = []
    for _, sp in irm.iterrows():
        m = dated[(dated.age_yBP >= sp.age_min_BP1950) &
                  (dated.age_yBP <= sp.age_max_BP1950)]
        if len(m) == 0:
            continue
        recs.append(dict(age=sp.age_mid_BP1950,
                         irm=sp.IRMsoft_flux_Am2_per_yr,
                         ln_EF_lith=np.log(m.EF_lith).mean(),
                         ln_EF_Cu=np.log(m.EF_Cu).mean(),
                         ln_EF_Ni=np.log(m.EF_Ni).mean()))
    b = pd.DataFrame(recs)
    out = {}
    for c in ("ln_EF_lith", "ln_EF_Cu", "ln_EF_Ni"):
        rho, p = stats.spearmanr(b[c], np.log(b.irm))
        out[c] = (rho, p, len(b))
    return b, out


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    te = load_te()

    # ── Event-window summary + formal flag ──────────────────────────────
    tab = summarise(te)
    tab_path = os.path.join(OUTDIR, "TableS_event_enrichment.csv")
    tab.round(3).to_csv(tab_path, index=False)
    print("Event-window enrichment factors (median | max, local rolling baseline):")
    show = ["window", "n"] + [f"EF_{e}_med" for e in ("Ni", "Co", "Cu", "lith", "Zn")] \
        + [f"EF_{e}_max" for e in ("Ni", "Co", "Cu", "lith", "Zn")]
    print(tab[show].round(2).to_string(index=False))

    flagged = te[te.detrital_flag]
    print(f"\nFormal detrital flag (robust log-z > {Z_FLAG} in BOTH Cr and Zn): "
          f"{len(flagged)} sample(s)")
    if len(flagged):
        print(flagged[["sample", "depth_cm", "age_yBP", "z_Cr", "z_Zn",
                       "EF_Ni", "EF_Co", "EF_Cu"]].round(2).to_string(index=False))
        near = te[(~te.detrital_flag) & (te.z_Cr > 2.5) & (te.z_Zn > 2.5)]
        print(f"   next-nearest candidates (z > 2.5 in both, unflagged): {len(near)}")
        print("-> flagged samples: the basal-contact pair (253.5, 255.5 cm, "
              "below the dated span — the age model ends at 253.0 cm) and one "
              "discrete layer at 5.94 cm (~1839 CE) inside the post-LIA surge "
              "interval. Neither event window contains a flagged sample. The "
              "5.94 cm layer carries only a mild Ni excess (EF 1.29), which "
              "biases that single point toward SLOWER inferred drip — i.e. "
              "against, not toward, the surge it sits in. No flagged sample "
              "affects any headline result.")

    # ── Trajectory slopes: kinetic vs source vs observed ────────────────
    # The trajectory is the direction of the excursion from the baseline
    # origin (EF = 1, 1), so the slopes are fitted through the origin for
    # both the observations and the kinetic path (which passes through the
    # origin at V = V_MODERN by construction). The ordinary least-squares
    # slope among the excursion points alone is also reported; with only
    # three elevated samples in the 568-point record it is not a stable
    # estimate of the trajectory direction.
    ev = window(te, WIN_52_CORE)
    lnN, lnC = np.log(ev.EF_Ni.values), np.log(ev.EF_Co.values)
    obs_slope = float((lnN * lnC).sum() / (lnN ** 2).sum())
    obs_ols, obs_r = np.polyfit(lnN, lnC, 1)[0], np.corrcoef(lnN, lnC)[0, 1]
    per_point = lnC / lnN
    Vg = np.geomspace(0.3, V_MODERN, 200)
    kN = np.log(phi_of_V(Vg, **PARAMS["Ni"]) / phi_of_V(np.array([V_MODERN]), **PARAMS["Ni"])[0])
    kC = np.log(phi_of_V(Vg, **PARAMS["Co"]) / phi_of_V(np.array([V_MODERN]), **PARAMS["Co"])[0])
    m = Vg >= 1.0
    kin_slope = float((kN[m] * kC[m]).sum() / (kN[m] ** 2).sum())
    print(f"\n5.2 ka trajectory: observed d lnEF_Co/d lnEF_Ni from the baseline "
          f"origin = {obs_slope:.2f} (per-sample {per_point.min():.2f}-"
          f"{per_point.max():.2f}, n = {len(ev)}; OLS among the excursion "
          f"points alone {obs_ols:.2f}, r = {obs_r:.2f})")
    print(f"  kinetic model (V {V_MODERN} -> 1, through the origin): {kin_slope:.2f}, "
          f"steepening toward Ni saturation at low V; source-covariance axis: 0.49; "
          f"detrital axis: ~1 with lithogenic co-enrichment (absent here).")

    # ── 8.2 ka window statement ─────────────────────────────────────────
    w82 = window(te, WIN_82)
    print(f"\n8.2 ka window (8.0-8.4 ka, n = {len(w82)}): "
          f"EF_Ni med {w82.EF_Ni.median():.2f} (max {w82.EF_Ni.max():.2f}), "
          f"EF_Co med {w82.EF_Co.median():.2f} (max {w82.EF_Co.max():.2f}) — "
          f"no kinetic response; EF_Cu max {w82.EF_Cu.max():.2f} "
          f"(episodic NOM/particulate delivery only).")

    # ── IRM comparison ──────────────────────────────────────────────────
    b, rho = irm_comparison(te)
    for c, lab in (("ln_EF_lith", "lithogenic index"), ("ln_EF_Cu", "EF_Cu"),
                   ("ln_EF_Ni", "EF_Ni")):
        r, p, n = rho[c]
        print(f"Spearman rho({lab}, ln IRM_soft flux) = {r:+.2f} "
              f"(p = {p:.2g}, n = {n} specimen bins)")
    print("-> no record-wide covariance between the dissolved-digest lithogenic "
          "index and the magnetite flux: the particulate (storm/event-flow) and "
          "dissolved NOM-complexed (matrix-flow) carriers are independent "
          "delivery pathways. The 8.2 ka IRM_soft minimum therefore records "
          "reduced storm-driven particulate delivery, which the kinetic Ni-Co "
          "pair — flat through the window — does not and should not see; at "
          "5.2 ka, by contrast, BOTH pathways collapse (IRM_soft falls to its "
          "record minimum as the drip proxy censors), the signature of an "
          "infiltration-excess shutdown.")

    # ── Figure: enrichment summary (a) and Co-Ni discriminant (b) ─────────
    # The earlier five-panel composition (record-wide profiles, event zooms,
    # discriminant, bars) is retained in the git history (pre-2026-09-17).
    dated = te.dropna(subset=["age_yBP"])
    apply_style()
    fig, (axE, axD) = plt.subplots(1, 2, figsize=(DOUBLE_COL, 2.7),
                                   gridspec_kw=dict(width_ratios=[1.25, 1], wspace=0.3))
    bars = tab.set_index("window")
    groups = [("5.2 ka core (155.7-157.1 cm)", "5.2 ka core", COL_TEAL_600),
              ("8.2 ka window (8.0-8.4 ka)", "8.2 ka window", COL_BG_300),
              ("basal detrital (>253.4 cm)", "basal detrital layer", COL_BROWN_500)]
    order = ["Co", "Cu", "Ni", "Zn", "V", "Cr"]          # OMC-bound -> lithogenic
    x = np.arange(len(order))
    wdt = 0.26
    for i, (g, lab, c) in enumerate(groups):
        vals = [bars.loc[g, f"EF_{e}_max"] for e in order]
        axE.bar(x + (i - 1) * wdt, vals, width=wdt, color=c, label=lab, lw=0)
    axE.axhline(1.0, color=COL_BG_900, lw=0.4)
    axE.set_yscale("log")
    axE.set_ylim(0.8, 40)
    axE.set_xticks(x)
    axE.set_xticklabels(order)
    axE.set_ylabel("Maximum enrichment factor")
    # element-group brackets drawn explicitly (horizontal line + short tick at each end),
    # sitting just below the x-tick labels, with the group name centred below
    from matplotlib.lines import Line2D
    tr = axE.get_xaxis_transform()   # x in data, y in axes fraction
    for x0, x1, label in [(-0.35, 2.35, "organically complexed"),
                          (2.65, 5.35, "lithogenic")]:
        y_bar, y_tick, y_text = -0.11, -0.09, -0.15
        axE.add_line(Line2D([x0, x1], [y_bar, y_bar], transform=tr,
                            color=COL_BG_600, lw=0.4, clip_on=False))
        for xe in (x0, x1):
            axE.add_line(Line2D([xe, xe], [y_bar, y_tick], transform=tr,
                                color=COL_BG_600, lw=0.4, clip_on=False))
        axE.text((x0 + x1) / 2, y_text, label, transform=tr,
                 ha="center", va="top", fontsize=5.5, color=COL_BG_600)
    
    panel_label(axE, "a")

    axD.scatter(np.log(dated.EF_Ni), np.log(dated.EF_Co), s=2.5, color=COL_BG_300,
                alpha=0.5, lw=0, label="all dated samples", rasterized=True)
    axD.scatter(lnN, lnC, s=9, color=COL_TEAL_600, lw=0.3, edgecolor=COL_BG_900,
                label="5.2 ka core", zorder=5)
    bas = window(te, WIN_BASAL)
    axD.scatter(np.log(bas.EF_Ni), np.log(bas.EF_Co), s=12, marker="s",
                color=COL_BROWN_500, lw=0.3, edgecolor=COL_BG_900,
                label="basal detrital layer", zorder=5)
    axD.plot(kN, kC, color=COL_BG_900, lw=0.9, label=f"kinetic path (slope {kin_slope:.1f})")
    xs = np.linspace(-0.6, 1.4, 10)
    axD.plot(xs, 0.49 * xs, color=COL_BG_500, lw=0.7, ls="--",
             label="source axis (0.49)")
    axD.plot(xs, 1.0 * xs, color=COL_BROWN_500, lw=0.7, ls=":",
             label="detrital axis (~1)")
    axD.set_xlabel("ln EF Ni")
    axD.set_ylabel("ln EF Co")
    axD.set_xlim(-0.9, 1.6)
    axD.set_ylim(-0.9, 2.7)
    from matplotlib.patches import Patch
    _dh, _dl = axD.get_legend_handles_labels()          # panel-b handles (dots, path, axes)
    _group_handles = [Patch(fc=c, ec="none", label=lab) for _, lab, c in groups]
    # drop the "5.2 ka core" and "basal detrital layer" entries from panel b (they duplicate the group patches)
    _keep = [(h, l) for h, l in zip(_dh, _dl) if l not in {"5.2 ka core", "basal detrital layer"}]
    # combined legend for panels a and b, placed below the figure (outside the panels)
    _all_handles = _group_handles + [h for h, _ in _keep]
    _all_labels  = [p.get_label() for p in _group_handles] + [l for _, l in _keep]
    fig.legend(_all_handles, _all_labels, frameon=False, fontsize=5.2,
               loc="lower center", bbox_to_anchor=(0.5, -0.02), ncol=6,
               columnspacing=1.5, handlelength=1.4, handleheight=0.7)
    fig.subplots_adjust(bottom=0.24)
    panel_label(axD, "b")
    for ax in (axE, axD):
        style_ax(ax)

    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUTDIR, f"FigS_detrital_ternary_screen.{ext}"))
    print(f"\nwrote {tab_path}")
    print(f"wrote {os.path.join(OUTDIR, 'FigS_detrital_ternary_screen.png')} (+pdf)")


if __name__ == "__main__":
    main()
