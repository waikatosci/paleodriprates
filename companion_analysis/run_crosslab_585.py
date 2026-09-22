#!/usr/bin/env python3
"""
run_crosslab_585.py -- builds the 585-point sensitivity input for Supplementary
Methods 14.3 / Supplementary Figure 17 and runs the native-resolution inversion
on it headlessly through Dr Paleo, at the canonical production settings
(dr_app/PRODUCTION_SETTINGS.md; parameter log in
manuscript_figures/external/run_logs/input_summary_native.csv).

Input 585 = the 568 primary-run rows (dr_app/HS4_example_inputs/HS4_TE_canonical.csv)
          + the 17 second-laboratory rows that are not depth duplicates of a
            primary-run row, rescaled with the paired-depth fit
            (primary = a * second + b, fitted at 8.09, 156.38 and 157.48 cm).

Usage
  python companion_analysis/run_crosslab_585.py            # 585-point run (sensitivity)
  python companion_analysis/run_crosslab_585.py --check    # 568-point run, compared with the released summary

Output
  manuscript_figures/external/drip_rate_summary_hr_corr585.csv   (--check writes nothing to external/)
  manuscript_figures/external/run_logs/input_summary_corr585.csv
"""
import argparse
import os
import shutil
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..'))
APP = os.path.join(ROOT, 'dr_app')
EXT = os.path.join(ROOT, 'manuscript_figures', 'external')
RAW = os.path.join(APP, 'HS4_example_inputs', 'HS4_TE.csv')
CANON = os.path.join(APP, 'HS4_example_inputs', 'HS4_TE_canonical.csv')
EXC = os.path.join(ROOT, 'calibration', 'excluded_points.csv')
PAIR_DEPTHS = [8.09, 156.38, 157.48]


def production_params(log_name: str = 'input_summary_native.csv'):
    """The canonical run's parameters, read from its own log (native by default;
    input_summary_1cm.csv and input_summary_ageprop.csv give the other two runs)."""
    log = pd.read_csv(os.path.join(EXT, 'run_logs', log_name))
    p = dict(zip(log.parameter, log.value.astype(str)))
    for k in ('run_id', 'timestamp'):
        p.pop(k, None)
    p['hiatus_zones'] = []
    p['temp_C'] = p.pop('cave_temperature_C')     # logged name -> run parameter name
    p['te_list'] = [{'col_depth': 'Depth (cm)', 'col_proxy': 'Ni (ppm)', 'unit': 'ppm'},
                    {'col_depth': 'Depth (cm)', 'col_proxy': 'Co (ppm)', 'unit': 'ppm'}]
    p['te1_col_depth'] = p['te2_col_depth'] = 'Depth (cm)'
    p['te1_col_proxy'], p['te2_col_proxy'] = 'Ni (ppm)', 'Co (ppm)'
    p['te1_unit'] = p['te2_unit'] = 'ppm'
    p['generate_realisations'] = False
    p['use_cached_proxy'] = False
    return p


def build_585():
    raw = pd.read_csv(RAW, encoding='utf-8-sig')
    d, ni, co = raw.columns[:3]
    exc = pd.read_csv(EXC)
    exc = exc[exc['ni_ppm'].notna()]
    second = np.zeros(len(raw), bool)
    for _, e in exc.iterrows():
        second |= (np.round(raw[d], 4) == round(float(e.depth_cm), 4)) & np.isclose(raw[ni], e.ni_ppm, atol=0.01)
    S = raw[second].copy()
    canon = pd.read_csv(CANON)
    canon.columns = [d, ni, co]
    fit = {}
    for col in (ni, co):
        xs, ys = [], []
        for dep in PAIR_DEPTHS:
            xs.append(S[np.isclose(S[d], dep, atol=0.006)].iloc[0][col])
            ys.append(canon[np.isclose(canon[d], dep, atol=0.006)].iloc[0][col])
        fit[col] = np.polyfit(xs, ys, 1)
        print(f'{col}: primary = {fit[col][0]:.4f} x second {fit[col][1]:+.4f}')
    keep = ~S[d].apply(lambda v: any(abs(v - p) < 0.006 for p in PAIR_DEPTHS))
    R = S[keep].copy()
    for col in (ni, co):
        R[col] = fit[col][0] * R[col] + fit[col][1]
    out = pd.concat([canon, R[[d, ni, co]]]).sort_values(d).reset_index(drop=True)
    print(f'input: {len(canon)} primary + {len(R)} rescaled second-laboratory rows = {len(out)}')
    return out


def run(te: pd.DataFrame, params: dict, workdir: str, extra_files: dict = None) -> pd.DataFrame:
    sys.path.insert(0, APP)
    sys.path.insert(0, ROOT)
    import app as drp
    up = os.path.join(workdir, 'uploads'); outd = os.path.join(workdir, 'outputs')
    for f in (up, outd):
        shutil.rmtree(f, ignore_errors=True); os.makedirs(f)
    te.to_csv(os.path.join(up, 'trace_elem1.csv'), index=False)
    for name, src in (extra_files or {}).items():
        shutil.copy(src, os.path.join(up, name))
    drp.UPLOAD_FOLDER, drp.OUTPUT_FOLDER = up, outd
    t0 = time.time()
    drp._run_model(dict(params))
    if drp.run_state.get('error'):
        raise SystemExit(drp.run_state['error'])
    print(f'run complete in {time.time() - t0:.0f} s')
    return pd.read_csv(os.path.join(outd, 'drip_rate_summary.csv')), outd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true', help='run the 568-point input and compare with the released summary')
    ap.add_argument('--workdir', default=os.path.join(ROOT, '.crosslab_run'))
    a = ap.parse_args()
    params = production_params()
    if a.check:
        te = pd.read_csv(CANON)
        s, _ = run(te, params, a.workdir)
        ref = pd.read_csv(os.path.join(EXT, 'drip_rate_summary_hr.csv'))
        m = s.merge(ref, on='depth', suffixes=('', '_ref'))
        dev = (m.pc50 / m.pc50_ref - 1).abs()
        print(f'568-point check: {len(m)}/{len(ref)} depths matched; max |pc50 deviation| = {100 * dev.max():.3f}%')
        return
    te = build_585()
    s, outd = run(te, params, a.workdir)
    s.to_csv(os.path.join(EXT, 'drip_rate_summary_hr_corr585.csv'), index=False)
    shutil.copy(os.path.join(outd, 'input_summary.csv'), os.path.join(EXT, 'run_logs', 'input_summary_corr585.csv'))
    ref = lambda df: df[((df.depth >= 140) & (df.depth < 155)) | ((df.depth > 158) & (df.depth <= 175))].pc50.median()
    w = s[(s.depth >= 150) & (s.depth <= 165)]; mn = w.loc[w.pc50.idxmin()]
    print(f'585-point run: {len(s)} depths; 5.2 ka minimum {mn.pc50:.2f} drips/min at {mn.depth:.2f} cm '
          f'({100 * (1 - mn.pc50 / ref(s)):.0f}% below the surrounding mid-Holocene median {ref(s):.1f})')


if __name__ == '__main__':
    main()
