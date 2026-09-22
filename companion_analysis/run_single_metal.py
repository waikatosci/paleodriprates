#!/usr/bin/env python3
"""
run_single_metal.py -- runs the native-resolution inversion on one proxy metal
at a time (Ni alone, Co alone), at the canonical production settings, so that
the drip rate each metal implies can be compared with the joint Ni-Co result
(Supplementary Methods 14.4; Supplementary Figure 19).

Input is the released 568-sample primary-run table
(dr_app/HS4_example_inputs/HS4_TE_canonical.csv); the run parameters are read
from the canonical run's own log (run_crosslab_585.production_params).

Usage
  python companion_analysis/run_single_metal.py --metal Ni [--mode native|1cm|ageprop|joint]
  python companion_analysis/run_single_metal.py --metal Co

--mode selects which canonical run is repeated (native resolution by default;
1cm = 1-cm-smoothed; ageprop = age-propagated on the U-Th age model).
--metal both reruns the joint two-metal inversion in that mode, as a check
against the released summary.

Output
  manuscript_figures/external/drip_rate_summary_<hr|lr|ap>_<metal>only.csv
  manuscript_figures/external/run_logs/input_summary_<mode>_<metal>only.csv
"""
import argparse
import os
import shutil

import pandas as pd

from run_crosslab_585 import APP, CANON, EXT, ROOT, production_params, run

LOGS = {'native': 'input_summary_native.csv', '1cm': 'input_summary_1cm.csv',
        'ageprop': 'input_summary_ageprop.csv'}
TAG = {'native': 'hr', '1cm': 'lr', 'ageprop': 'ap'}
AGE_DEPTH = os.path.join(APP, 'HS4_example_inputs', 'HS4_age_depth.csv')


def single_metal_params(metal: str, mode: str = 'native') -> dict:
    p = production_params(LOGS[mode])
    if mode == 'ageprop':
        p['col_depth'] = 'Distance (cm)'
        p['col_age'] = 'Age(corr) (yrs BP)'
        p['col_age_err'] = 'Age error (yrs BP 2-sigma)'
    if metal == 'both':
        return p
    src = 'te1_' if metal == 'Ni' else 'te2_'
    for key in list(p):
        if key.startswith(src):
            p['te1_' + key[len(src):]] = p[key]
    for key in [k for k in p if k.startswith('te2_')]:
        p.pop(key)
    col = f'{metal} (ppm)'
    p['te_list'] = [{'col_depth': 'Depth (cm)', 'col_proxy': col, 'unit': 'ppm'}]
    p['te1_col_depth'], p['te1_col_proxy'], p['te1_unit'] = 'Depth (cm)', col, 'ppm'
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--metal', choices=['Ni', 'Co', 'both'], required=True)
    ap.add_argument('--mode', choices=['native', '1cm', 'ageprop'], default='native')
    ap.add_argument('--workdir', default=None)
    a = ap.parse_args()
    workdir = a.workdir or os.path.join(ROOT, f'.single_{a.metal}_{a.mode}_run')
    te = pd.read_csv(CANON)
    te.columns = ['Depth (cm)', 'Ni (ppm)', 'Co (ppm)']
    if a.metal != 'both':
        te = te[['Depth (cm)', f'{a.metal} (ppm)']]
    extra = {'depth_age.csv': AGE_DEPTH} if a.mode == 'ageprop' else None
    s, outd = run(te, single_metal_params(a.metal, a.mode), workdir, extra)
    x = 'age' if a.mode == 'ageprop' else 'depth'
    if a.metal == 'both':
        ref = pd.read_csv(os.path.join(EXT, f'drip_rate_summary_{TAG[a.mode]}.csv'))
        m = s.merge(ref, on=x, suffixes=('', '_ref'))
        print(f'joint {a.mode} check: {len(m)}/{len(ref)} rows matched; max |pc50 deviation| = '
              f'{100 * (m.pc50 / m.pc50_ref - 1).abs().max():.3f}%')
        return
    tag = f'{a.metal}only'
    s.to_csv(os.path.join(EXT, f'drip_rate_summary_{TAG[a.mode]}_{tag}.csv'), index=False)
    shutil.copy(os.path.join(outd, 'input_summary.csv'),
                os.path.join(EXT, 'run_logs', f'input_summary_{a.mode}_{tag}.csv'))
    if a.mode == 'ageprop':
        # 5.2 ka window on the age axis; surrounding mid-Holocene background 4.5-6.0 ka outside 5.0-5.4 ka
        ev = s[(s.age >= 5000) & (s.age <= 5400)]
        bg = s[((s.age >= 4500) & (s.age < 5000)) | ((s.age > 5400) & (s.age <= 6000))]
        mn = ev.loc[ev.pc50.idxmin()]
        print(f'{a.metal}-only {a.mode}: 5.2 ka minimum {mn.pc50:.2f} drips/min at {mn.age:.0f} yr BP; '
              f'{100 * (1 - mn.pc50 / bg.pc50.median()):.0f}% below the surrounding median {bg.pc50.median():.1f}')
        return
    ref = s[((s.depth >= 140) & (s.depth < 155)) | ((s.depth > 158) & (s.depth <= 175))].pc50.median()
    w = s[(s.depth >= 150) & (s.depth <= 165)]
    mn = w.loc[w.pc50.idxmin()]
    print(f'{a.metal}-only {a.mode}: {len(s)} depths; 5.2 ka minimum {mn.pc50:.2f} drips/min at {mn.depth:.2f} cm '
          f'(IQR {mn.pc25:.2f}-{mn.pc75:.2f}); {100 * (1 - mn.pc50 / ref):.0f}% below the surrounding '
          f'mid-Holocene median {ref:.1f}')
    print(w[['depth', 'pc25', 'pc50', 'pc75']].round(2).to_string(index=False))


if __name__ == '__main__':
    main()
