#!/usr/bin/env python3
"""
crosslab_replication.py -- the second-laboratory run (November 2019) as an independent
replication of the 5.2 ka event on Ni (Supplementary Methods 14.3; Supplementary Figure 19b).

Each 2019 value is corrected in two steps:
  1. for the Ca of its own digest, using the Ca slope fitted across the 13 samples outside
     the 5.2 ka core (crosslab_matrix_check.ca_matrix_baseline), relative to the run mean;
  2. onto the primary-run scale, with the linear relation between the two runs at the three
     depths both laboratories analysed (8.09, 156.38, 157.48 cm).
The 17 corrected samples that are not depth duplicates are added to the 568 primary-run samples
and the Ni-only inversion is run at the production settings (run_single_metal.py). Ni is used
because Co carries the surplus described in Supplementary Methods 14.4.

Output
  manuscript_figures/external/drip_rate_summary_hr_Nionly_rep585.csv
  manuscript_figures/external/run_logs/input_summary_native_Nionly_rep585.csv
  manuscript_figures/external/crosslab_replication_2019_corrected.csv
"""
import contextlib
import io
import os
import shutil

import numpy as np
import pandas as pd

import crosslab_matrix_check as cm
import run_crosslab_585 as rc
from run_single_metal import single_metal_params

PAIRS = {8.09: 'HS4-C-987', 156.38: 'HS4-A-885', 157.48: 'HS4-A-875'}
BG = lambda d: ((d >= 140) & (d < 155)) | ((d > 158) & (d <= 175))


def corrected_2019():
    lab2 = cm.load_lab2(cm.DEFAULT_IN)
    with contextlib.redirect_stdout(io.StringIO()):
        base = cm.ca_matrix_baseline(lab2)
    canon = pd.read_csv(rc.CANON)
    canon.columns = ['depth', 'Ni', 'Co']
    out = lab2[['depth_cm', 'Ca', 'Ni', 'Co']].copy()
    for el in ('Ni', 'Co'):
        v = lab2[el] - base[el]['slope_per_ppm_Ca'] * (lab2['Ca'] - lab2['Ca'].mean())
        x = [v[s] for s in PAIRS.values()]
        y = [canon.loc[np.isclose(canon.depth, d, atol=0.006), el].iloc[0] for d in PAIRS]
        a, c = np.polyfit(x, y, 1)
        print(f'{el}: Ca-normalised second-run value -> primary scale: {a:.3f} x {c:+.3f} '
              f'(r = {np.corrcoef(x, y)[0, 1]:.3f} at the three shared depths)')
        out[el + '_corrected'] = a * v + c
    return out, canon


def main():
    corr, canon = corrected_2019()
    corr.to_csv(os.path.join(rc.EXT, 'crosslab_replication_2019_corrected.csv'))
    extra = corr[~corr.index.isin(PAIRS.values())]
    te = pd.concat([canon[['depth', 'Ni']],
                    pd.DataFrame({'depth': extra.depth_cm, 'Ni': extra.Ni_corrected})]).sort_values('depth')
    te.columns = ['Depth (cm)', 'Ni (ppm)']
    with contextlib.redirect_stdout(io.StringIO()):
        s, outd = rc.run(te, single_metal_params('Ni'), os.path.join(rc.ROOT, '.replication_run'))
    s.to_csv(os.path.join(rc.EXT, 'drip_rate_summary_hr_Nionly_rep585.csv'), index=False)
    shutil.copy(os.path.join(outd, 'input_summary.csv'),
                os.path.join(rc.EXT, 'run_logs', 'input_summary_native_Nionly_rep585.csv'))
    shutil.rmtree(os.path.join(rc.ROOT, '.replication_run'), ignore_errors=True)
    ref = s.loc[BG(s.depth), 'pc50'].median()
    w = s[(s.depth >= 154.5) & (s.depth <= 158.4)].copy()
    w['source'] = np.where(w.depth.round(2).isin(extra.depth_cm.round(2)), 'second run', 'primary run')
    mn = w.loc[w.pc50.idxmin()]
    prim = pd.read_csv(os.path.join(rc.EXT, 'drip_rate_summary_hr_Nionly.csv'))
    pm = prim[(prim.depth >= 154.5) & (prim.depth <= 158.4)].pc50.min()
    print(f'Ni-only, 568 + 17 corrected second-run samples: minimum {mn.pc50:.2f} drips/min at {mn.depth:.2f} cm '
          f'({mn.source}), {100 * (1 - mn.pc50 / ref):.0f}% below the surrounding median {ref:.1f} '
          f'(primary run alone: {pm:.2f})')
    print(w[['depth', 'source', 'pc25', 'pc50', 'pc75']].round(2).to_string(index=False))
    ns = s[s.depth <= 12]
    print(f'near-surface check (0-12 cm): minimum {ns.pc50.min():.1f}, median {ns.pc50.median():.1f} drips/min')


if __name__ == '__main__':
    main()
