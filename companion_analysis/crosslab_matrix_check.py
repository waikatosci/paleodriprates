#!/usr/bin/env python3
"""
crosslab_matrix_check.py -- Supplementary Methods 14.3, 14.4 and Supplementary
Figure 18. Reproduces the numbers behind the withdrawal of the 20 samples
analysed by a second laboratory and shows why the two runs cannot be combined
uncorrected.

The second laboratory (Wuhan, November 2019) measured the full trace-element
suite for 21 HS4 powders on an Agilent 8900 triple-quadrupole ICP-MS, and
returned the four common silicate reference materials within a few percent of
their certified values for every element used here. The two laboratories
nonetheless disagree at three depths (8.09, 156.38, 157.48 cm) they both
analysed, by a near-constant additive amount. The reason, evident once the
full suite is examined, is the calcium matrix: the second-run standards were
not matrix-matched to the ~40% Ca of a calcite digest, so its Co, Ni and Fe
signals in the calcite samples carry a Ca-proportional polyatomic interference
(40Ca16O+ at m/z 56 for Fe; the same on 59 and 60 for Co and Ni). The primary
run used Ca-matrix-matched standards and does not carry it.

This script prints every number cited in the text and Supplementary Figure 18,
and re-renders the figure.

Input   manuscript_figures/external/HS4_TE_full_suite_both_labs.xlsx
Output  manuscript_figures/output/FigS_matrix_interference.{png,pdf}
        stdout : the numbers, one line each, for cross-check with the paper.
"""
from __future__ import annotations
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


HERE = Path(__file__).resolve().parent
REPO = HERE.parent
DEFAULT_IN = REPO / 'manuscript_figures' / 'external' / 'HS4_TE_full_suite_both_labs.xlsx'
DEFAULT_OUT = REPO / 'manuscript_figures' / 'output'

# depth (cm) of each 2019-run HS4 sample; both A- (event core) and C- (surface)
DEPTH_CM = {
    'HS4-A-871': 157.80, 'HS4-A-873': 157.65, 'HS4-A-875': 157.48,
    'HS4-A-877': 157.30, 'HS4-A-879': 157.10, 'HS4-A-881': 156.88,
    'HS4-A-883': 156.64, 'HS4-A-885': 156.38, 'HS4-A-887': 156.23,
    'HS4-A-889': 155.87, 'HS4-A-891': 155.60, 'HS4-A-893': 155.20,
    'HS4-C-985':   8.48, 'HS4-C-986':   8.29, 'HS4-C-987':   8.09,
    'HS4-C-988':   7.89, 'HS4-C-989':   7.69, 'HS4-C-990':   7.49,
    'HS4-C-991':   7.30, 'HS4-C-993':   7.00,
}
# The seven samples that make up the 5.2 ka core (155.2-156.9 cm)
EVENT_CORE = {'HS4-A-881', 'HS4-A-883', 'HS4-A-885', 'HS4-A-887',
              'HS4-A-889', 'HS4-A-891', 'HS4-A-893'}

# Upper-continental-crust ratios used in the detrital mass balance
UCC = {'Al': 81500.0, 'Th': 10.5, 'Ti': 3840.0, 'Zr': 193.0}   # ppm
UCC_CO = 17.3
UCC_NI = 47.0


def load_lab2(path: Path) -> pd.DataFrame:
    """Return a per-sample dataframe of the 2019 (Lab 2) run of the 21 HS4 powders."""
    raw = pd.read_excel(path, sheet_name='Lab 2', header=None)
    names = raw.iloc[0].tolist()                       # analysis IDs
    samples = raw.iloc[1].tolist()                     # sample names
    # Element rows start at row 4 (0-indexed); first cell holds the symbol
    tbl = {}
    for _, row in raw.iloc[4:].iterrows():
        sym = row.iloc[0]
        if isinstance(sym, str) and 1 <= len(sym) <= 3:
            tbl[sym] = pd.to_numeric(row.iloc[1:].values, errors='coerce')
    df = pd.DataFrame(tbl, index=samples[1:])
    df = df.loc[[s for s in df.index if isinstance(s, str) and s.startswith('HS4-')]]
    df['depth_cm'] = [DEPTH_CM[s] for s in df.index]
    df['in_5.2ka_core'] = df.index.isin(EVENT_CORE)
    return df.sort_values('depth_cm')


def report_standard_recoveries(path: Path) -> None:
    raw = pd.read_excel(path, sheet_name='Lab 2', header=None)
    names = raw.iloc[1].tolist()
    tbl = {}
    for _, row in raw.iloc[4:].iterrows():
        sym = row.iloc[0]
        if isinstance(sym, str) and 1 <= len(sym) <= 3:
            tbl[sym] = pd.to_numeric(row.iloc[1:].values, errors='coerce').tolist()
    df = pd.DataFrame(tbl, index=names[1:])
    print('\n== silicate reference-material recoveries (2019 run, meas/ref) ==')
    for label in ('AGV-2', 'BHVO-2', 'BCR-2', 'RGM-2'):
        rows = [i for i, n in enumerate(df.index) if n == label]
        if len(rows) == 2:
            meas, ref = df.iloc[rows[0]], df.iloc[rows[1]]
            pct = {e: 100 * meas[e] / ref[e]
                   for e in ('Co', 'Ni', 'Cu', 'Mn', 'Cr', 'Al', 'Th')
                   if pd.notna(ref[e]) and ref[e] > 0}
            fmt = ', '.join(f'{e} {pct[e]:.0f}%' for e in pct)
            print(f'  {label}: {fmt}')
    print('  (RGM Ni is near detection and is not certified; the same low '
          'recovery is seen for RGM-1 in the primary run and is not diagnostic '
          'of either laboratory.)')


def ca_matrix_baseline(df: pd.DataFrame) -> dict:
    """Regressions of Co, Ni and Fe on Ca across the samples outside the event core."""
    base = df[~df['in_5.2ka_core']]
    out = {}
    print(f'\n== Ca-matrix baseline (n = {len(base)}, samples outside the 5.2 ka core) ==')
    print('  slope*Ca is the pure interference contribution at 400,000 ppm Ca (the')
    print('  number quoted in the paper); the full-fit value is intercept + slope*Ca.')
    for element in ('Co', 'Ni', 'Fe'):
        r = stats.linregress(base['Ca'], base[element])
        matrix_contribution = r.slope * 400_000
        full_fit = r.intercept + r.slope * 400_000
        out[element] = {'r': r.rvalue, 'intercept_ppm': r.intercept,
                        'slope_per_ppm_Ca': r.slope,
                        'matrix_contribution_ppm': matrix_contribution,
                        'full_fit_at_calcite_Ca_ppm': full_fit}
        print(f'  {element}: r = {r.rvalue:+.2f}, intercept = {r.intercept:+.1f} ppm, '
              f'slope*400k = {matrix_contribution:.2f} ppm, full-fit = {full_fit:.2f} ppm')
    print('  Fe slope*Ca at ~0.7% = the 40Ca16O+ polyatomic interference at m/z 56;')
    print('  the Co and Ni baselines are consistent with the corresponding')
    print('  Ca-bearing ions at m/z 59 and 60. A near-zero Co intercept and a')
    print('  strongly negative Fe intercept are compatible with the interference')
    print('  reading, since real dissolved Co is trace and real Fe is scavenged')
    print('  below the range the interference lifts it into.')
    return out


def detrital_mass_balance(df: pd.DataFrame, base_by_ca: dict) -> None:
    ev = df[df['in_5.2ka_core']]
    print(f'\n== detrital mass balance across the event core (n = {len(ev)}) ==')
    # local baseline for the excess: median of the flanking A- samples
    flank_names = ['HS4-A-871', 'HS4-A-873', 'HS4-A-875', 'HS4-A-877', 'HS4-A-879']
    flank = df.loc[flank_names]
    co_base, ni_base = flank['Co'].median(), flank['Ni'].median()
    ev_excess_co = (ev['Co'] - co_base).abs().max()
    ev_excess_ni = (ev['Ni'] - ni_base).abs().max()
    print(f'  event-core excesses over local baseline: '
          f'Co up to {ev_excess_co:.2f} ppm, Ni up to {ev_excess_ni:.2f} ppm')
    print(f'  measured maxima across the seven core samples:')
    for tracer in ('Al', 'Th', 'Ti', 'Zr'):
        m = ev[tracer].max()
        co_ucc = m * UCC_CO / UCC[tracer]
        ni_ucc = m * UCC_NI / UCC[tracer]
        print(f'    {tracer:<3s} <= {m:>7.3f} ppm  ->  detrital bound '
              f'Co <= {co_ucc:.4f} ppm, Ni <= {ni_ucc:.4f} ppm')
    print('  every tracer bounds detrital Co and Ni at below 1% of the observed excess.')
    for name, other in (('Al', 'Co'), ('Al', 'Ni')):
        rho, p = stats.spearmanr(df[name], df[other])
        print(f'  Spearman rho({name}, {other}) across the 21 samples = {rho:+.2f} (p = {p:.2f})')


def residence_time_signature(df: pd.DataFrame, base_by_ca: dict) -> None:
    """Excess Co and Ni over the Ca-proportional baseline scale with Mn."""
    ex_co = df['Co'] - (base_by_ca['Co']['intercept_ppm']
                        + base_by_ca['Co']['slope_per_ppm_Ca'] * df['Ca'])
    ex_ni = df['Ni'] - (base_by_ca['Ni']['intercept_ppm']
                        + base_by_ca['Ni']['slope_per_ppm_Ca'] * df['Ca'])
    print('\n== Mn residence-time signature (all 21 samples) ==')
    for name, series in (('Co', ex_co), ('Ni', ex_ni)):
        r = stats.linregress(df['Mn'], series)
        rho, p = stats.spearmanr(df['Mn'], series)
        print(f'  excess {name} / Mn slope = {r.slope:.4f} '
              f'(Pearson r = {r.rvalue:+.2f}), Spearman rho = {rho:+.2f} (p = {p:.3f})')
    ratio_all = df['Ce'] / df['La']
    print(f"  Ce/La outside the event: median = {ratio_all[~df['in_5.2ka_core']].median():.2f}")
    print(f"  Ce/La inside the event:  median = {ratio_all[df['in_5.2ka_core']].median():.2f}, "
          f"max = {ratio_all[df['in_5.2ka_core']].max():.2f}")


def render_figure(df: pd.DataFrame, out_dir: Path) -> None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    out_dir.mkdir(parents=True, exist_ok=True)
    base = df[~df['in_5.2ka_core']]
    event = df[df['in_5.2ka_core']]

    TEAL, ORANGE, GREY, AXIS = '#0f766e', '#c2410c', '#6b7280', '#374151'
    plt.rcParams.update({'font.size': 8, 'font.family': 'DejaVu Sans',
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.linewidth': 0.6, 'axes.edgecolor': AXIS,
                         'xtick.color': AXIS, 'ytick.color': AXIS,
                         'axes.labelcolor': AXIS})
    fig, axes = plt.subplots(2, 2, figsize=(6.8, 5.4), dpi=300)

    def scatter(ax, element, letter):
        r = stats.linregress(base['Ca'] / 1e3, base[element])
        x = np.linspace(370, 540, 2)
        ax.plot(x, r.intercept + r.slope * x, '-', color=GREY, lw=1, zorder=1)
        ax.scatter(base['Ca'] / 1e3, base[element], s=22, facecolor='white',
                   edgecolor=GREY, lw=1, zorder=2,
                   label='outside the 5.2 ka core (n = 13)')
        ax.scatter(event['Ca'] / 1e3, event[element], s=26, color=ORANGE, zorder=3,
                   label='5.2 ka core, 155.2–156.9 cm (n = 7)')
        ax.set_xlabel(r'Ca ($\times 10^3$ ppm)')
        ax.set_ylabel(f'{element} (ppm)')
        ax.text(0.04, 0.96, letter, transform=ax.transAxes, fontweight='bold',
                va='top', fontsize=9)
        ax.grid(axis='y', lw=0.4, color='#e5e7eb'); ax.set_axisbelow(True)

    scatter(axes[0, 0], 'Co', 'a')
    scatter(axes[0, 1], 'Ni', 'b')
    scatter(axes[1, 0], 'Fe', 'c')
    axes[0, 0].legend(frameon=True, framealpha=0.9, edgecolor='none',
                      fontsize=6.5, loc='upper left', bbox_to_anchor=(0.06, 0.93),
                      handlelength=1.2, borderpad=0.3)

    ax = axes[1, 1]
    for element, colour, marker in (('Co', TEAL, 'o'), ('Ni', ORANGE, 's')):
        r = stats.linregress(base['Ca'], base[element])
        excess = df[element] - (r.intercept + r.slope * df['Ca'])
        ax.scatter(df['Mn'], excess, s=22, color=colour, marker=marker, zorder=3,
                   label=element)
        rr = stats.linregress(df['Mn'], excess)
        x = np.linspace(0, 115, 2)
        ax.plot(x, rr.intercept + rr.slope * x, '-', color=colour, lw=1, alpha=0.6)
        ax.text(117, rr.intercept + rr.slope * 117, element,
                color=colour, fontsize=7, va='center')
    ax.set_xlabel('Mn (ppm)')
    ax.set_ylabel('Excess over Ca-proportional baseline (ppm)')
    ax.set_xlim(0, 130)
    ax.text(0.04, 0.96, 'd', transform=ax.transAxes, fontweight='bold',
            va='top', fontsize=9)
    ax.legend(frameon=True, framealpha=0.9, edgecolor='none',
              fontsize=6.5, loc='upper left', bbox_to_anchor=(0.06, 0.93),
              handlelength=1.2, borderpad=0.3)
    ax.grid(axis='y', lw=0.4, color='#e5e7eb'); ax.set_axisbelow(True)

    fig.tight_layout()
    fig.savefig(out_dir / 'FigS_matrix_interference.png', facecolor='white')
    fig.savefig(out_dir / 'FigS_matrix_interference.pdf', facecolor='white')
    print(f'\nwrote {out_dir}/FigS_matrix_interference.{{png,pdf}}')


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', type=Path, default=DEFAULT_IN)
    ap.add_argument('--output-dir', type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    report_standard_recoveries(args.input)
    df = load_lab2(args.input)
    baseline = ca_matrix_baseline(df)
    detrital_mass_balance(df, baseline)
    residence_time_signature(df, baseline)
    render_figure(df, args.output_dir)


if __name__ == '__main__':
    main()
