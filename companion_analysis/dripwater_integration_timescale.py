#!/usr/bin/env python3
"""
dripwater_integration_timescale.py -- over what timescale does the metal supply at the HS4 drip
track drip rate? (Supplementary Methods 16.)

The kinetic proxy is attributed to baseflow, the median flow of a hydrological year, while calcite
is deposited preferentially in summer. This script asks the monitoring record how much of the
seasonal cycle survives in the water that feeds the stalagmite: it compares the seasonal amplitude
of drip rate with that of the metal/Ca ratios, and correlates the ratios with the mean drip rate
over windows from one day to two years.

Input
  dripwater/HS4_dripwater_canonical.csv   HS4 drip monitoring, 2007-2015 (date, drip rate, solution
                                          chemistry; Me/Ca in umol/mol)

Output
  manuscript_figures/output/TableS_dripwater_integration.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / 'manuscript_figures' / 'output'
WINDOWS = ['1D', '30D', '90D', '180D', '365D', '730D']

d = (pd.read_csv(REPO / 'dripwater' / 'HS4_dripwater_canonical.csv', parse_dates=['date'])
       .rename(columns={'DR_drops_per_min': 'DR', 'NiCa_umol_per_mol': 'NiCa',
                        'CoCa_umol_per_mol': 'CoCa', 'Ca_ppb': 'Ca'})
       .dropna(subset=['DR']).sort_values('date').set_index('date'))
print(f'{len(d)} samples, {d.index.min().date()} to {d.index.max().date()}, '
      f'{d.index.year.nunique()} calendar years')

# seasonal amplitude: range of the monthly means, as a fraction of the overall mean
rows = []
for c in ['DR', 'NiCa', 'CoCa', 'Ca']:
    m = d.groupby(d.index.month)[c].mean()
    rows.append(dict(quantity=c, n=int(d[c].notna().sum()),
                     seasonal_amplitude=round((m.max() - m.min()) / m.mean(), 3),
                     cv=round(d[c].std() / d[c].mean(), 3)))
seasonal = pd.DataFrame(rows)
print('\nSeasonal amplitude (range of monthly means / mean):')
print(seasonal.to_string(index=False))

# how well does the mean drip rate over a window explain the metal ratios?
rows = []
for w in WINDOWS:
    mean_dr = d['DR'].rolling(w).mean()
    row = dict(window=w)
    for c in ['NiCa', 'CoCa']:
        k = mean_dr.notna() & d[c].notna()
        r, p = stats.pearsonr(mean_dr[k], d[c][k])
        rho, _ = stats.spearmanr(mean_dr[k], d[c][k])
        row[f'r_{c}'] = round(r, 2)
        row[f'p_{c}'] = float(f'{p:.1e}')
        row[f'rho_{c}'] = round(rho, 2)
        row['n'] = int(k.sum())
    rows.append(row)
tab = pd.DataFrame(rows)
print('\nCorrelation of dripwater Me/Ca with mean drip rate over each window:')
print(tab.to_string(index=False))
print('\nNOTE: the record spans eight years, so the long-window coefficients rest on few independent '
      'intervals. The monotonic rise with window length, not any single coefficient, is the result.')

s = d.groupby(d.index.month)['DR'].mean()
summer, winter = s.loc[6:9].mean(), s.loc[[12, 1, 2, 3]].mean()
print(f'\nDrip rate: summer (Jun-Sep) {summer:.1f}, winter (Dec-Mar) {winter:.1f} drips/min '
      f'(ratio {summer / winter:.2f}); annual baseflow mean used for the calibration 16.66, '
      f'summer mean {100 * (summer / 16.66 - 1):+.0f}% against it')

OUT.mkdir(exist_ok=True)
out = OUT / 'TableS_dripwater_integration.csv'
with open(out, 'w') as f:
    seasonal.to_csv(f, index=False)
    f.write('\n')
    tab.to_csv(f, index=False)
    f.write(f'\nsummer_mean_drips_min,{summer:.2f}\nwinter_mean_drips_min,{winter:.2f}\n'
            f'calibration_annual_baseflow_mean,16.66\n')
print(f'\nsaved {out}')
