#!/usr/bin/env python3
"""
element_covariance.py -- which trace elements covary with Ni and Co, in the calcite record,
in the full-suite runs at the 5.2 ka horizon, and in the monitored HS4 dripwater.

1. Primary-run record (manuscript_figures/external/HS4_TE_multielement.csv; Ni, Co, Cu, Cr, V, Zn;
   n = 568): Spearman rho on log concentrations, raw and after removing a 61-sample rolling
   median (the Supplementary Methods 12 baseline), with p-values for an effective sample size
   corrected for lag-1 autocorrelation; lagged cross-correlation of the detrended series
   (-5 to +5 samples); enrichment of each element at the three 5.2 ka samples.
2. 2019 run (HS4_TE_full_suite_both_labs.xlsx, sheet Lab 2): Spearman rho of every element with
   the Co and Ni excess over the Ca-proportional baseline of Supplementary Methods 14.3, across
   the 12 samples at 155.2-157.8 cm and within the 7 core samples.
3. Dripwater (dripwater/HS4_dripwater_canonical.csv; n = 92, 2007-2016): Spearman rho of every
   element and of drip rate with Ni and Co.

Output: manuscript_figures/output/element_covariance_*.csv and FigS_element_covariance.{png,pdf}
"""
import contextlib
import io
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
EXT = REPO / 'manuscript_figures' / 'external'
OUT = REPO / 'manuscript_figures' / 'output'
EVENT = [155.74, 156.38, 157.01]
REC_EL = ['Ni', 'Co', 'Cu', 'Cr', 'V', 'Zn']


def spearman_neff(a, b):
    """Spearman rho with a two-sided p-value for an effective n reduced by lag-1 autocorrelation."""
    m = a.notna() & b.notna()
    a, b = a[m].reset_index(drop=True), b[m].reset_index(drop=True)
    r = stats.spearmanr(a, b)[0]
    ra, rb = a.rank().autocorr(1), b.rank().autocorr(1)
    ne = len(a) * (1 - ra * rb) / (1 + ra * rb)
    t = r * np.sqrt((ne - 2) / (1 - r * r))
    return r, 2 * stats.t.sf(abs(t), ne - 2), len(a), ne


def record():
    d = pd.read_csv(EXT / 'HS4_TE_multielement.csv').dropna(subset=['Ni', 'Co'])
    d = d.sort_values('depth_cm').reset_index(drop=True)
    ev = pd.Series(np.isclose(d.depth_cm.values[:, None], EVENT, atol=0.01).any(1))
    L = np.log(d[REC_EL].where(d[REC_EL] > 0))
    Z = L - L.rolling(61, center=True, min_periods=15).median()
    rows, lag_rows = [], []
    for ref in ('Ni', 'Co'):
        for e in REC_EL:
            if e == ref:
                continue
            for name, X, keep in (('raw', L, None), ('detrended', Z, None),
                                  ('detrended, 5.2 ka removed', Z, ~ev)):
                a, b = (X[ref], X[e]) if keep is None else (X[ref][keep], X[e][keep])
                r, p, n, ne = spearman_neff(a, b)
                rows.append(dict(ref=ref, element=e, series=name, rho=r, p_adj=p, n=n, n_eff=ne))
            for k in range(-5, 6):
                b = Z[e].shift(k); m = Z[ref].notna() & b.notna() & ~ev
                lag_rows.append(dict(ref=ref, element=e, lag=k, rho=stats.spearmanr(Z[ref][m], b[m])[0]))
    ef = np.exp(Z[ev]).assign(depth_cm=d.depth_cm[ev].values)
    return pd.DataFrame(rows), pd.DataFrame(lag_rows), ef


def suite2019():
    sys.path.insert(0, str(HERE))
    import crosslab_matrix_check as cm
    df = cm.load_lab2(cm.DEFAULT_IN)
    with contextlib.redirect_stdout(io.StringIO()):
        base = cm.ca_matrix_baseline(df)
    for e in ('Co', 'Ni'):
        df['x' + e] = df[e] - (base[e]['intercept_ppm'] + base[e]['slope_per_ppm_Ca'] * df['Ca'])
    A = df[df.index.str.startswith('HS4-A')]
    core = A[A['in_5.2ka_core']]
    els = [c for c in df.columns if c not in ('depth_cm', 'in_5.2ka_core', 'xCo', 'xNi', 'Co', 'Ni')
           and df[c].notna().all() and df[c].std() > 0]
    rows = []
    for e in els:
        rows.append(dict(element=e,
                         rho_Co_12=stats.spearmanr(A.xCo, A[e])[0], rho_Ni_12=stats.spearmanr(A.xNi, A[e])[0],
                         rho_Co_core7=stats.spearmanr(core.xCo, core[e])[0], rho_Ni_core7=stats.spearmanr(core.xNi, core[e])[0],
                         EF_core_median=core[e].median() / A[~A['in_5.2ka_core']][e].median()))
    return pd.DataFrame(rows).set_index('element').sort_values('rho_Co_core7', ascending=False)


def dripwater():
    d = pd.read_csv(REPO / 'dripwater' / 'HS4_dripwater_canonical.csv')
    cols = [c for c in d.columns if c.endswith('_ppb')] + ['DR_drops_per_min']
    rows = []
    for c in cols:
        row = dict(element=c.replace('_ppb', '').replace('DR_drops_per_min', 'drip rate'))
        for ref in ('Ni_ppb', 'Co_ppb'):
            if c == ref:
                continue
            m = d[ref].notna() & d[c].notna() & (d[c] > 0)
            row['rho_' + ref[:2]], row['p_' + ref[:2]] = stats.spearmanr(d.loc[m, ref], d.loc[m, c])
            row['n'] = int(m.sum())
        rows.append(row)
    return pd.DataFrame(rows).set_index('element')


def figure(rec, lag, ef, drip):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    NI, CO, AXIS, GRID = '#1f5fa8', '#6b7280', '#374151', '#e5e7eb'
    plt.rcParams.update({'font.size': 8, 'font.family': 'DejaVu Sans', 'axes.spines.top': False,
                         'axes.spines.right': False, 'axes.linewidth': 0.6, 'axes.edgecolor': AXIS,
                         'xtick.color': AXIS, 'ytick.color': AXIS, 'axes.labelcolor': AXIS})
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.9), dpi=300, gridspec_kw={'width_ratios': [1, 1, 1.25]})
    # a: record-wide detrended rho
    a = ax[0]; order = ['Co', 'Ni', 'Cr', 'Cu', 'Zn', 'V']
    det = rec[rec.series == 'detrended']
    for i, e in enumerate(order):
        for ref, col, mk, off in (('Ni', NI, 'o', -0.12), ('Co', CO, 's', 0.12)):
            r = det[(det.ref == ref) & (det.element == e)]
            if len(r):
                a.plot(r.rho.iloc[0], i + off, marker=mk, color=col, ms=5, ls='none',
                       label=f'with {ref}' if i == (0 if ref == 'Ni' else 1) else None)
    a.axvline(0, color=AXIS, lw=0.6); a.set_yticks(range(len(order))); a.set_yticklabels(order); a.invert_yaxis()
    a.set_xlim(-0.2, 0.6); a.set_xlabel('Spearman ρ (detrended, n = 568)')
    a.grid(axis='x', lw=0.4, color=GRID); a.set_axisbelow(True)
    a.legend(frameon=False, fontsize=6.5, loc='lower right', handletextpad=0.3)
    a.text(-0.3, 1.03, 'a', transform=a.transAxes, fontweight='bold', fontsize=9)
    # b: enrichment at the three 5.2 ka samples
    b = ax[1]; els = ['Ni', 'Co', 'Cu', 'Cr', 'V', 'Zn']
    for i, e in enumerate(els):
        v = ef[e].values
        col = NI if e == 'Ni' else (CO if e == 'Co' else '#9ca3af')
        b.plot(v, [i] * len(v), 'o', color=col, ms=4.5, alpha=0.9)
        b.plot([v.min(), v.max()], [i, i], '-', color=col, lw=1)
    b.axvline(1, color=AXIS, lw=0.6); b.set_yticks(range(len(els))); b.set_yticklabels(els); b.invert_yaxis()
    b.set_xlabel('Enrichment at the 5.2 ka samples (×)'); b.set_xlim(0.5, 5)
    b.grid(axis='x', lw=0.4, color=GRID); b.set_axisbelow(True)
    b.text(-0.3, 1.03, 'b', transform=b.transAxes, fontweight='bold', fontsize=9)
    # c: dripwater rho with Ni and Co
    c = ax[2]; show = ['Co', 'Ni', 'Pb', 'La', 'Al', 'Zn', 'Cu', 'Mn', 'Cd', 'drip rate', 'Mg']
    for i, e in enumerate(show):
        for ref, col, mk, off in (('Ni', NI, 'o', -0.12), ('Co', CO, 's', 0.12)):
            key = 'rho_' + ref
            if e in drip.index and key in drip.columns and pd.notna(drip.loc[e, key]):
                c.plot(drip.loc[e, key], i + off, marker=mk, color=col, ms=5, ls='none')
    c.axvline(0, color=AXIS, lw=0.6); c.set_yticks(range(len(show))); c.set_yticklabels(show); c.invert_yaxis()
    c.set_xlim(-0.8, 0.9); c.set_xlabel('Spearman ρ in dripwater (n = 92)')
    c.grid(axis='x', lw=0.4, color=GRID); c.set_axisbelow(True)
    c.text(-0.3, 1.03, 'c', transform=c.transAxes, fontweight='bold', fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / 'FigS_element_covariance.png', facecolor='white')
    fig.savefig(OUT / 'FigS_element_covariance.pdf', facecolor='white')


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rec, lag, ef = record()
    s19 = suite2019()
    drip = dripwater()
    rec.to_csv(OUT / 'element_covariance_record.csv', index=False)
    lag.to_csv(OUT / 'element_covariance_record_lagged.csv', index=False)
    ef.to_csv(OUT / 'element_covariance_5p2ka_enrichment.csv', index=False)
    s19.to_csv(OUT / 'element_covariance_2019run.csv')
    drip.to_csv(OUT / 'element_covariance_dripwater.csv')
    pd.set_option('display.width', 200)
    print('== record, Spearman rho with Ni and Co ==')
    print(rec.pivot_table(index=['ref', 'element'], columns='series', values='rho').round(2).to_string())
    print('\n== detrended, p-values for effective n ==')
    print(rec[rec.series == 'detrended'][['ref', 'element', 'rho', 'p_adj', 'n_eff']].round(3).to_string(index=False))
    print('\n== lagged (detrended, 5.2 ka removed), lags -5..+5 ==')
    for (ref, e), g in lag.groupby(['ref', 'element']):
        print(f'  {ref}-{e}: ' + ' '.join(f'{v:+.2f}' for v in g.sort_values('lag').rho))
    print('\n== enrichment at the three 5.2 ka samples (value / 61-sample rolling median) ==')
    print(ef.round(2).to_string(index=False))
    print('\n== 2019 run: rho with Co and Ni excess (12 samples, and the 7 core samples) ==')
    print(s19.round(2).to_string())
    print('\n== dripwater: rho with Ni and Co ==')
    print(drip.sort_values('rho_Ni', ascending=False).round(3).to_string())
    figure(rec, lag, ef, drip)
    print(f'\nwrote {OUT}/FigS_element_covariance.{{png,pdf}} and element_covariance_*.csv')


if __name__ == '__main__':
    main()
