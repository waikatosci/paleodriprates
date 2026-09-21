#!/usr/bin/env python3
"""
events_differ_lithogenic.py -- reproduces every number in Supplementary Methods 13.2-13.3
(dust and PCP at 8.2 ka; their absence at 5.2 ka) from released data.

Inputs
  manuscript_figures/external/HS4_8p2ka_digest_2017.xlsx  sheet 'Processed'
      51 carbonate-digest ICP-MS samples across 230.8-236.8 cm (8.06-8.46 ka), full element
      suite (Al, K, Mg, Sr, Ca, Cr, Co, Ni ...). First column block = digest ppb; ratios use Ca 43.
  manuscript_figures/external/HS4_TE_multielement.csv      whole-record Cr, V (589 samples)
  manuscript_figures/HS4_SourceData.xlsx                   03_isotopes + 02_chronology

Definitions
  8.2 ka interval  = 230.8-236.8 cm (the digest-suite depth window)
  5.2 ka core      = 155.2-157.8 cm (as in Supplementary Methods 12)
  whole-record contrasts: median in window / median of all other samples, Mann-Whitney U
  isotope contrasts: window mean minus mean of a +/-0.5 ka flanking baseline, Mann-Whitney U
  detrital mass balance: all Al assumed detrital with upper-continental-crust composition
      (Rudnick & Gao 2003: Al 8.15 wt%, Co 17.3 ppm, Ni 47 ppm)
"""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
EXT = os.path.join(ROOT, 'manuscript_figures', 'external'); WB = os.path.join(ROOT, 'manuscript_figures', 'HS4_SourceData.xlsx')
W82 = (230.8, 236.8); W52 = (155.2, 157.8)

DIGEST = os.path.join(EXT, 'HS4_8p2ka_digest_2017.xlsx')
if os.path.exists(DIGEST):
    d = pd.read_excel(DIGEST, sheet_name='Processed', header=None)
    hdr = [str(h) for h in d.iloc[3]]; b = d.iloc[4:].copy(); b.columns = hdr
    b = b.loc[:, ~b.columns.duplicated()]; b = b[b['Depth'].notna()]
    for c in b.columns:
        if c not in ('Sample ID', 'Aliquot', 'ID'): b[c] = pd.to_numeric(b[c], errors='coerce')
else:
    # the same 51-sample digest block is released as workbook sheet 07_lithogenic_8p2ka
    b = pd.read_excel(WB, sheet_name='07_lithogenic_8p2ka', header=3)
    b = b.rename(columns={'depth_cm': 'Depth', 'age_yBP': 'Age', 'Ca_ppb': 'Ca 43', 'Mg_ppb': 'Mg 24',
                          'Al_ppb': 'Al 27', 'K_ppb': 'K 39', 'Cr_ppb': 'Cr 52', 'Co_ppb': 'Co 59',
                          'Ni_ppb': 'Ni 60', 'Sr_ppb': 'Sr 88'})
    b = b[b['Depth'].notna()]
    for c in b.columns:
        if c != 'sample': b[c] = pd.to_numeric(b[c], errors='coerce')
    print('(digest read from HS4_SourceData.xlsx sheet 07_lithogenic_8p2ka)')
b['MgCa'] = b['Mg 24'] / b['Ca 43']; b['SrCa'] = b['Sr 88'] / b['Ca 43']; b['AlCa'] = b['Al 27'] / b['Ca 43']
print(f"8.2 ka digest suite: n = {len(b)}, {b.Depth.min():.1f}-{b.Depth.max():.1f} cm, {b.Age.min():.0f}-{b.Age.max():.0f} yr BP")

def rho(a, c):
    x = b[[a, c]].dropna(); r, p = spearmanr(x[a], x[c]); return r, p, len(x)
print("\nS13.2 Al covariance (Spearman):")
for a, c, lab in [('Al 27', 'K 39', 'Al-K'), ('Al 27', 'Cr 52', 'Al-Cr'), ('Al 27', 'Co 59', 'Al-Co'), ('Al 27', 'Ni 60', 'Al-Ni')]:
    r, p, n = rho(a, c); print(f"  rho({lab}) = {r:+.2f}  p = {p:.2g}  n = {n}")
detCo = (b['Al 27'] * 17.3 / 81500 / b['Co 59']).median(); detNi = (b['Al 27'] * 47 / 81500 / b['Ni 60']).median()
print(f"  detrital mass balance (UCC): Co {100*detCo:.2f}%  Ni {100*detNi:.2f}%  (median across window)")

print("\nS13.3 PCP:")
for a, c, lab in [('MgCa', 'SrCa', 'Mg/Ca-Sr/Ca'), ('AlCa', 'MgCa', 'Al/Ca-Mg/Ca'), ('AlCa', 'SrCa', 'Al/Ca-Sr/Ca')]:
    r, p, n = rho(a, c); print(f"  rho({lab}) = {r:+.2f}  p = {p:.2g}  n = {n}")

te = pd.read_csv(os.path.join(EXT, 'HS4_TE_multielement.csv'))
def contrast(col, lo, hi):
    m = (te.depth_cm >= lo) & (te.depth_cm <= hi); a = te[m][col].dropna(); o = te[~m][col].dropna()
    return a.median() / o.median(), mannwhitneyu(a, o)[1], len(a)
print("\nWhole-record contrasts (window median / rest-of-record median, Mann-Whitney):")
for col in ('Cr', 'V'):
    for lab, w in (('8.2 ka', W82), ('5.2 ka', W52)):
        f, p, n = contrast(col, *w); print(f"  {col} {lab}: x{f:.2f}  p = {p:.2g}  n = {n}")

iso = pd.read_excel(WB, sheet_name='03_isotopes', header=4).apply(pd.to_numeric, errors='coerce').dropna(subset=['depth_cm'])
ch = pd.read_excel(WB, sheet_name='02_chronology', header=6)[['depth_cm', 'age_yBP']].apply(pd.to_numeric, errors='coerce').dropna().sort_values('depth_cm')
iso['age'] = np.interp(iso.depth_cm, ch.depth_cm, ch.age_yBP)
print("\nIsotope contrasts (window mean - +/-0.5 ka flanking mean, Mann-Whitney):")
for lab, (lo, hi) in (('8.2 ka', W82), ('5.2 ka', W52)):
    ev = iso[(iso.depth_cm >= lo) & (iso.depth_cm <= hi)]; alo, ahi = ev.age.min(), ev.age.max()
    fl = iso[((iso.age >= alo - 500) & (iso.age < alo)) | ((iso.age > ahi) & (iso.age <= ahi + 500))]
    for col, name in (('d13C_permil', 'd13C'), ('d18O_permil', 'd18O')):
        print(f"  {lab} {name}: {ev[col].mean()-fl[col].mean():+.2f} permil  p = {mannwhitneyu(ev[col], fl[col])[1]:.2g}  (n = {len(ev)} vs {len(fl)})")
