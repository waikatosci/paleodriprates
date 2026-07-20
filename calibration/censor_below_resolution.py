"""
censor_below_resolution.py — flag Dr Paleo drip points that fell below the
joint Ni-Co resolution limit and emit a censored series.

WHY THIS EXISTS (2026-07-20)
Under the canonical kinetic-window dissociation width (sigma = pi/sqrt(6)),
the joint drip posterior is the geometric mean of the single-proxy Ni and Co
PDFs. Where the two proxies indicate strongly discordant rates — which happens
only at the extreme slow-drip minimum of the record — the geometric-mean
product collapses and Dr Paleo returns an all-zero column (every percentile 0).
This is NOT a bug and NOT a settings error: the old broad sigma = 1.39 kernel
smeared the two proxies together and manufactured a finite (wrong) value there;
the sharper, correct kernel honestly reports "below resolution."

For HS4 Run 1 (native canonical) this affects exactly two points, 156.23 and
156.64 cm, inside the megadrought minimum. Their per-proxy implied rates are
Ni ~2.8 and Co ~0.4 drips/min; the geometric mean (~1.0) sits just below the
coupled inversion's resolution limit. The adjacent survivor at 156.88 cm
resolves at 1.07 drips/min but with a boundary-collapsed posterior (relative
width ~1.7% vs the record-typical ~20%) — i.e. it too is pinned at the floor.
The empirical resolution limit is therefore ~1 drip/min; the sensitivity floor
(where the forward curve flattens against the equilibrium ceiling) is lower
still (~0.18 drips/min), so this is a joint-consistency limit, not saturation.

WHAT IT DOES
- Detects censored points: rows where ALL percentile columns are 0.
- Emits <stem>_censored.csv with the raw percentiles blanked to NaN on censored
  rows, plus three added columns:
    censored           1 on censored rows, 0 otherwise
    drip_upper_bound   = FLOOR on censored rows (the reportable "<" bound), NaN otherwise
    pc50_for_continuity= FLOOR on censored rows, pc50 otherwise
                         (use THIS for plotting/RQA/spectral so a literal 0 does
                          not create a spurious spike; FLOOR imputation is
                          conservative — it makes the trough look no deeper than
                          the limit, never deeper than the data support)
- Prints a run summary and a ready-to-paste Methods sentence.
- With --plot, renders a diagnostic of the affected interval (resolved points
  with 5-95% bars; censored points as open markers at the floor with a "<"; the
  floor line). This is a WORKING diagnostic, not an ngeo-style final figure.

FLOOR
Default 1.0 drip/min — the clean reportable limit. The empirical lowest-
resolved value is 1.07 (documented above); 1.0 is the round, slightly
conservative bound. Override with --floor if a run's floor differs.

USAGE
    python3 censor_below_resolution.py --input <run>/drip_rate_summary.csv \
        --floor 1.0 --out drip_rate_summary_censored.csv --plot censor_diag.png
"""
import argparse
import os

import numpy as np
import pandas as pd

PC_COLS = ["pc05", "pc10", "pc25", "pc50", "pc75", "pc90", "pc95"]


def censor(df, floor):
    pcs = [c for c in PC_COLS if c in df.columns]
    is_cens = (df[pcs] == 0).all(axis=1)
    out = df.copy()
    out["censored"] = is_cens.astype(int)
    out["drip_upper_bound"] = np.where(is_cens, floor, np.nan)
    out["pc50_for_continuity"] = np.where(is_cens, floor, out["pc50"])
    out.loc[is_cens, pcs] = np.nan          # blank raw percentiles on censored rows
    return out, df.loc[is_cens]


def methods_sentence(cens_df, floor):
    depths = ", ".join(f"{d:.2f}" for d in cens_df["depth"].tolist())
    return (
        "Drip rates were quantified from the joint posterior of the Ni and Co "
        "dissociation models (geometric mean of the two single-proxy PDFs). "
        "Under the kinetic-window dissociation width (\u03c3 = \u03c0/\u221a6), "
        "this joint posterior collapses where the two proxies indicate strongly "
        "discordant rates, which occurs only at the extreme minimum of the "
        f"record. At {len(cens_df)} points within the drought minimum "
        f"({depths} cm depth), inferred drip rates fell below the joint "
        f"resolution limit of \u2248{floor:g} drip min\u207b\u00b9 (the lowest "
        "rate the coupled inversion resolves) and are reported as censored "
        f"(\u2264{floor:g} drip min\u207b\u00b9); they are retained in the "
        "drip-rate series as upper bounds and excluded from spectral and "
        "recurrence analyses."
    )


def diagnostic_plot(full, cens, floor, path, pad_cm=8.0):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    lo = cens["depth"].min() - pad_cm
    hi = cens["depth"].max() + pad_cm
    win = full[(full["depth"] >= lo) & (full["depth"] <= hi)]
    res = win[win["pc50"] > 0]
    cen = win[win["pc50"] == 0]

    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.axhline(floor, ls="--", lw=1, color="#c0392b",
               label=f"resolution limit \u2248{floor:g} drip min\u207b\u00b9")
    yerr = np.vstack([res["pc50"] - res["pc05"], res["pc95"] - res["pc50"]])
    ax.errorbar(res["depth"], res["pc50"], yerr=yerr, fmt="o", ms=4,
                color="#2c3e50", ecolor="#95a5a6", elinewidth=0.8, capsize=2,
                label="resolved (5\u201395%)")
    ax.scatter(cen["depth"], np.full(len(cen), floor), s=70,
               facecolors="none", edgecolors="#c0392b", linewidths=1.6, zorder=5,
               label="censored (\u2264 limit)")
    for d in cen["depth"]:
        ax.annotate("", xy=(d, floor * 0.55), xytext=(d, floor),
                    arrowprops=dict(arrowstyle="->", color="#c0392b", lw=1.4))
        ax.text(d, floor * 1.06, "<", ha="center", va="bottom",
                color="#c0392b", fontsize=11, fontweight="bold")
    ax.set_xlabel("Depth (cm)")
    ax.set_ylabel("Drip rate (drips min\u207b\u00b9)")
    ax.set_title("Drought minimum: resolved vs censored drip rates "
                 "(working diagnostic)")
    ax.set_ylim(0, max(res["pc95"].max() * 1.05, floor * 2))
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax.grid(alpha=0.15)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Dr Paleo *_summary.csv")
    ap.add_argument("--floor", type=float, default=1.0,
                    help="resolution floor in drips/min (default 1.0)")
    ap.add_argument("--out", default=None, help="output censored CSV path")
    ap.add_argument("--plot", default=None, help="optional diagnostic PNG path")
    args = ap.parse_args()

    df = pd.read_csv(args.input)
    out, cens = censor(df, args.floor)

    n = len(cens)
    print(f"input: {args.input}  ({len(df)} rows)")
    print(f"floor: {args.floor:g} drips/min")
    print(f"censored points: {n}")
    if n:
        for _, r in cens.iterrows():
            print(f"  depth {r['depth']:.2f} cm  (was all-zero \u2192 \u2264{args.floor:g})")

    out_path = args.out or os.path.splitext(args.input)[0] + "_censored.csv"
    out.to_csv(out_path, index=False)
    print(f"wrote {out_path}")

    if n:
        print("\n--- Methods sentence ---")
        print(methods_sentence(cens, args.floor))

    if args.plot and n:
        diagnostic_plot(df, cens, args.floor, args.plot)
        print(f"\nwrote diagnostic {args.plot}")


if __name__ == "__main__":
    main()
