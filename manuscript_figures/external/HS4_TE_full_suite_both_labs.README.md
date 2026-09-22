# HS4_TE_full_suite_both_labs.xlsx

Solid-phase ICP-MS analyses of HS4 stalagmite powders measured for the full
trace-element suite (Li to U): 20 powders in the 2019 second-laboratory run
and five in the 2009 primary run. Circulated by C. Hu on 21 September 2026
as the primary evidence for the decisions taken in Supplementary Methods
14.3 and 14.4. Provided as Source Data with the paper, with C. Hu's permission
(22 September 2026).

## Sheets

**Lab 2** — Second-laboratory supplementary run, November 2019 (Wuhan
SampleSolution Analytical Technology Co., Ltd; Agilent 7700e ICP-MS, same
method as the 2009 run; analysis prefix W19111...).
20 HS4 samples: 12 from the 5.2 ka interval (HS4-A-871…893, 155.2–157.8 cm)
and 8 from near the growth surface (HS4-C-985…993, 7.0–8.5 cm), plus a
procedural blank (HS1-0), four silicate reference materials (AGV-2, BHVO-2,
BCR-2, RGM-2) and two laboratory blanks.

**lab 1** — Primary run, 2009 (Agilent 7500a; the run used throughout the
paper), full suite for five powders across the event (A-870/875/880/885/890; two of
the three depths shared with the 2019 run, 156.38 and 157.48 cm, are among
them, and A-880 is the resolved 5.2 ka minimum at 157.01 cm), and standards.

Concentrations are ppm in the calcite (silicate reference materials also in
ppm). Row `<u:</u> 10-9 (ng/g)…` (last row of sheet Lab 2) records the
distinction for the blanks and standards.

## What the two runs show

* Both are correctly calibrated against silicate reference materials
  (AGV-2 / BHVO-2 / BCR-2 return Co at 98–101% and Ni at 95–101% of the
  reference values in the 2019 run; comparable in the 2009 run).
* Across the 13 second-run HS4 samples that lie outside the 5.2 ka core,
  Co, Ni and Fe are linear functions of the Ca in the digest (Pearson
  r = 0.98, 0.96, 0.97; near-zero intercepts; ~1.6 ppm Co, ~5.8 ppm Ni,
  ~6,900 ppm Fe at calcite Ca). Fe at ~0.7% is the ⁴⁰Ca¹⁶O⁺ polyatomic
  interference at m/z 56; the Co and Ni baselines are consistent with the
  corresponding Ca-bearing ions at m/z 59 and 60. The primary run's
  standards were matrix-matched to high-purity Ca and do not carry it.
* Across the seven core samples, Al ≤38 ppm, Th ≤0.013 ppm, Ti ≤7.1 ppm and
  Zr ≤0.13 ppm. At upper-continental-crust ratios this bounds detrital
  Co at ≤0.03 ppm and Ni at ≤0.09 ppm: sample by sample, 1.0–3.5% of the
  excess Co and 1.6–6.0% of the excess Ni, and 1–2% for the samples that
  carry the largest excesses. Neither Co nor Ni correlates with Al across
  the 20 samples (Spearman ρ = 0.16 and 0.27).
* The correlations with Ca are anchored by the sample with the highest
  reported Ca (HS4-A-873, 536 ×10³ ppm); without it they remain
  significant (r = 0.88 Co, 0.68 Ni, 0.80 Fe; p ≤ 0.014, n = 12), the Co
  and Fe slopes are unchanged and the Ni slope falls to ~4.2 ppm at
  calcite Ca. Reported Ca runs 376–536 ×10³ ppm and exceeds the
  stoichiometric 400 ×10³ ppm of calcite in nine samples, so the Ca axis is
  an analytical quantity of each digest, not a solid composition.
* Once the Ca-proportional baseline is removed, the excess Co and Ni at
  the 5.2 ka horizon scale with Mn (Co/Mn = 0.022, Ni/Mn = 0.039;
  Spearman ρ(Co,Mn) = 0.85, ρ(Ni,Mn) = 0.70), Mn rises from ~1 ppm to
  25–111 ppm and Cu from ~0.6 to 1.5–2.9 ppm. The light rare earths are
  enriched roughly tenfold in the core (Ce 0.14–0.65 ppm against ≤0.05 ppm
  outside it), consistent with a minor Mn-oxide component, but Ce/La in the
  core (1.0–2.3) overlaps its range outside the core (0.8–1.7) and does not
  track Mn, so a cerium anomaly is not resolved.

## Reproduced by

`companion_analysis/crosslab_matrix_check.py` prints all numbers cited in
Supplementary Methods 14.3, 14.4 and Supplementary Figure 18, and re-renders
Supplementary Figure 18 from this file.
