#!/usr/bin/env python3
"""
ni_co_agreement.py -- compares the drip rate implied by Ni alone with that implied
by Co alone, sample by sample, and reports the 5.2 ka interval
(Supplementary Methods 14.4; Supplementary Figure 19).

Inputs (all in manuscript_figures/external/):
  drip_rate_summary_hr.csv          joint Ni-Co inversion (released record)
  drip_rate_summary_hr_Nionly.csv   Ni alone   } run_single_metal.py --metal Ni
  drip_rate_summary_hr_Coonly.csv   Co alone   } run_single_metal.py --metal Co
  drip_rate_summary_hr_Nionly_rep585.csv and crosslab_replication_2019_corrected.csv
                                    the corrected second-run samples (crosslab_replication.py)
  ../../dr_app/HS4_example_inputs/HS4_TE_canonical.csv   the 568 primary-run samples

The carrier bounds solve, for each 5.2 ka sample, for the drip rate at which the
kinetic model accounts for the measured Ni and Co once a Co-bearing phase with a
fixed Ni/Co mass ratio has been removed. The forward model is the one inverted by
Dr Paleo (labile fraction 1 - (1 - lambda_F) E[exp(-k tau)], tau = 60/d, log-normal
k with the production mu and sigma), anchored so that the record-median
concentration maps to the calibration drip rate; it reproduces the single-metal
posterior medians to within 0.5%.

Output
  manuscript_figures/output/FigS_ni_co_agreement.{png,pdf}
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq

REPO = Path(__file__).resolve().parent.parent
EXT = REPO / 'manuscript_figures' / 'external'
OUT = REPO / 'manuscript_figures' / 'output'
CANON = REPO / 'dr_app' / 'HS4_example_inputs' / 'HS4_TE_canonical.csv'
EVENT = [155.74, 156.38, 157.01]
BG = lambda d: ((d >= 140) & (d < 155)) | ((d > 158) & (d <= 175))

SIGMA, LAMBDA_F, D_CAL = np.pi / np.sqrt(6), 0.01, 16.66
MU = {'Ni': -3.8372, 'Co': -5.4171}
# record-median calcite concentrations the calibration anchors to D_CAL (Supplementary Table 2)
MEDIAN = {'Ni': 3.6382, 'Co': 0.3121}
_k = np.linspace(-20, 20, 8001)


def load():
    j = pd.read_csv(EXT / 'drip_rate_summary_hr.csv')
    ni = pd.read_csv(EXT / 'drip_rate_summary_hr_Nionly.csv')
    co = pd.read_csv(EXT / 'drip_rate_summary_hr_Coonly.csv')
    te = pd.read_csv(CANON)
    te.columns = ['depth', 'Ni', 'Co']
    m = (te.merge(j[['depth', 'pc25', 'pc50', 'pc75']], on='depth')
           .merge(ni[['depth', 'pc25', 'pc50', 'pc75']].add_suffix('_ni').rename(columns={'depth_ni': 'depth'}), on='depth')
           .merge(co[['depth', 'pc25', 'pc50', 'pc75']].add_suffix('_co').rename(columns={'depth_co': 'depth'}), on='depth'))
    return m


def labile(metal, d):
    w = np.exp(-(_k - MU[metal]) ** 2 / (2 * SIGMA ** 2)); w /= w.sum()
    return 1 - (1 - LAMBDA_F) * (w * np.exp(-np.exp(_k) * 60.0 / d)).sum()


def forward(metal, d, median):
    return median * labile(metal, d) / labile(metal, D_CAL)


def report(m):
    ev = m.depth.isin(EVENT)
    lr = np.log(m.pc50_ni / m.pc50_co)
    q05, q50, q95 = np.exp(lr[~ev].quantile([.05, .5, .95]))
    z = (lr[ev] - lr[~ev].mean()) / lr[~ev].std()
    print(f'Ni-alone / Co-alone drip rate, the {int((~ev).sum())} samples outside 5.2 ka: '
          f'median {q50:.2f}, 5th-95th percentile {q05:.2f}-{q95:.2f}')
    print('5.2 ka samples: ' + '; '.join(
        f'{r.depth:.2f} cm Ni {r.pc50_ni:.1f}, Co {r.pc50_co:.1f}, joint {r.pc50:.1f} '
        f'(ratio {r.pc50_ni / r.pc50_co:.1f}, z = {zz:.1f}, rank {int(lr.rank(ascending=False)[i])})'
        for (i, r), zz in zip(m[ev].iterrows(), z)))
    for col, name in (('pc50', 'joint'), ('pc50_ni', 'Ni alone'), ('pc50_co', 'Co alone')):
        ref = m.loc[BG(m.depth), col].median()
        mn = m.loc[ev, col].min()
        print(f'  {name:9s}: minimum {mn:.2f} drips/min, {100 * (1 - mn / ref):.0f}% below the surrounding median {ref:.1f}')
    faster = m[['pc50_ni', 'pc50_co']].max(axis=1)
    lo = m.loc[faster.nsmallest(3).index]
    print('lowest drip rates on which both metals agree (the faster of the two):',
          ', '.join(f'{r.depth:.2f} cm {f:.1f}' for (_, r), f in zip(lo.iterrows(), faster[lo.index])))
    med = MEDIAN
    check = [forward('Ni', r.pc50_ni, med['Ni']) / r.Ni - 1 for _, r in m[ev].iterrows()]
    print(f'forward-model check at the 5.2 ka samples: max |Ni mismatch| {100 * max(map(abs, check)):.1f}%')
    for rc in (0.0, 1.0, 1.8):
        parts = []
        for _, r in m[ev].iterrows():
            f = lambda d: (r.Ni - forward('Ni', d, med['Ni'])) - rc * (r.Co - forward('Co', d, med['Co']))
            d = brentq(f, 0.3, 200)
            parts.append(f'{r.depth:.2f} cm {d:.1f} (Co-bearing phase {100 * (1 - forward("Co", d, med["Co"]) / r.Co):.0f}% of Co)')
        print(f'  carrier Ni/Co = {rc}: ' + '; '.join(parts))


def figure(m):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    NI, CO, ORANGE, GREY, AXIS = '#1f5fa8', '#6b7280', '#c2410c', '#9ca3af', '#374151'
    plt.rcParams.update({'font.size': 8, 'font.family': 'DejaVu Sans',
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.linewidth': 0.6, 'axes.edgecolor': AXIS,
                         'xtick.color': AXIS, 'ytick.color': AXIS, 'axes.labelcolor': AXIS})
    ev = m.depth.isin(EVENT)
    q05, q95 = np.exp(np.log(m.pc50_ni / m.pc50_co)[~ev].quantile([.05, .95]))
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.3), dpi=300, gridspec_kw={'width_ratios': [1, 1.35]})
    x = np.array([1, 100])
    a.fill_between(x, x * q05, x * q95, color=GREY, alpha=0.2, lw=0, zorder=0,
                   label='5th–95th percentile elsewhere')
    a.plot(x, x, '-', color=AXIS, lw=0.8, zorder=1)
    a.scatter(m.pc50_co[~ev], m.pc50_ni[~ev], s=9, facecolor='white', edgecolor=GREY, lw=0.6, zorder=2,
              label=f'other samples (n = {int((~ev).sum())})')
    a.scatter(m.pc50_co[ev], m.pc50_ni[ev], s=26, color=ORANGE, zorder=3, label='5.2 ka samples (n = 3)')
    for _, r in m[ev].iterrows():
        a.annotate(f'{r.depth:.2f} cm', (r.pc50_co, r.pc50_ni), xytext=(-6, -2), ha='right',
                   textcoords='offset points', fontsize=6.5, color=ORANGE)
    from matplotlib.ticker import FixedLocator, NullFormatter, ScalarFormatter
    a.set_xscale('log'); a.set_yscale('log'); a.set_xlim(1.0, 90); a.set_ylim(1.5, 90)
    for axis in (a.xaxis, a.yaxis):
        axis.set_major_locator(FixedLocator([1, 2, 5, 10, 20, 50]))
        axis.set_major_formatter(ScalarFormatter()); axis.set_minor_formatter(NullFormatter())
    a.set_xlabel(r'Drip rate from Co alone (drips min$^{-1}$)')
    a.set_ylabel(r'Drip rate from Ni alone (drips min$^{-1}$)')
    a.text(50, 58, '1:1', fontsize=6.5, ha='right', color=AXIS)
    a.legend(frameon=False, fontsize=6.3, loc='upper left', handlelength=1.2)
    a.text(-0.22, 1.02, 'a', transform=a.transAxes, fontweight='bold', fontsize=9)

    w = m[(m.depth >= 148) & (m.depth <= 165)]
    b.axvspan(155.5, 157.2, color=ORANGE, alpha=0.07, lw=0)
    for sfx, col, mk, lab in (('_ni', NI, 'o', 'Ni alone'), ('_co', CO, 's', 'Co alone')):
        b.fill_between(w.depth, w['pc25' + sfx], w['pc75' + sfx], color=col, alpha=0.15, lw=0)
        b.plot(w.depth, w['pc50' + sfx], '-', marker=mk, color=col, ms=2.8, lw=1, label=lab)
    b.plot(w.depth, w.pc50, '--', color='#111827', lw=1, label='joint Ni–Co (published record)')
    rep = EXT / 'drip_rate_summary_hr_Nionly_rep585.csv'
    corr = EXT / 'crosslab_replication_2019_corrected.csv'
    if rep.exists() and corr.exists():
        r5 = pd.read_csv(rep); c19 = pd.read_csv(corr)
        dep = c19.depth_cm[~c19.depth_cm.round(2).isin([8.09, 156.38, 157.48])]
        r5 = r5[r5.depth.round(2).isin(dep.round(2)) & (r5.depth >= 148) & (r5.depth <= 165)]
        b.errorbar(r5.depth, r5.pc50, yerr=[r5.pc50 - r5.pc25, r5.pc75 - r5.pc50], fmt='D', ms=3.6,
                   mfc='white', mec=ORANGE, ecolor=ORANGE, elinewidth=0.8, capsize=0, zorder=4,
                   label='Ni, second run (corrected)')
    b.set_xlim(165, 148); b.set_ylim(0, 40)
    b.set_xlabel('Depth (cm)'); b.set_ylabel(r'Drip rate (drips min$^{-1}$)')
    b.legend(frameon=False, fontsize=6.3, loc='upper right', ncol=2, handlelength=1.6, columnspacing=1.0)
    r = m[m.depth == 157.01].iloc[0]
    b.annotate(f'157.01 cm: Ni {r.pc50_ni:.1f}, Co {r.pc50_co:.1f},\njoint {r.pc50:.1f}',
               (157.01, r.pc50_ni), xytext=(164.6, 2.2), fontsize=6.5,
               arrowprops=dict(arrowstyle='-', lw=0.5, color=AXIS))
    b.text(-0.12, 1.02, 'b', transform=b.transAxes, fontweight='bold', fontsize=9)
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / 'FigS_ni_co_agreement.png', facecolor='white')
    fig.savefig(OUT / 'FigS_ni_co_agreement.pdf', facecolor='white')
    print(f'wrote {OUT}/FigS_ni_co_agreement.{{png,pdf}}')


if __name__ == '__main__':
    m = load()
    report(m)
    figure(m)
