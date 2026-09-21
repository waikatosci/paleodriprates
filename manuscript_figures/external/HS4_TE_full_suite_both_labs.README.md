# HS4_TE_full_suite_both_labs.xlsx

Solid-phase ICP-MS analyses of 21 HS4 stalagmite powders, measured for the
full trace-element suite (Li to U) by both analytical runs used in this
paper. Circulated by C. Hu on 21 September 2026 as the primary evidence for
the decisions taken in Supplementary Methods 14.3 and 14.4.

## Sheets

**Lab 2** — Second-laboratory supplementary run, November 2019 (Agilent
8900 triple-quadrupole ICP-MS, analysis prefix W19111...).
21 HS4 samples: 12 from the 5.2 ka interval (HS4-A-871…893, 155.2–157.8 cm)
and 8 from near the growth surface (HS4-C-985…993, 7.0–8.5 cm), plus four
silicate reference materials (AGV-2, BHVO-2, BCR-2, RGM-2) and two blanks.

**lab 1** — Primary run, 2009 (Agilent 7500a; the run used throughout the
paper), same three shared depths only: A-870/875/880/885/890, and standards.

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
* Across the 21 samples, Al ≤38 ppm, Th ≤0.013 ppm, Ti ≤7.1 ppm and
  Zr ≤0.13 ppm. At upper-continental-crust ratios this bounds detrital
  Co at ≤0.03 ppm and Ni at ≤0.09 ppm — below 1% of the measured event
  excesses.
* Once the Ca-proportional baseline is removed, the excess Co and Ni at
  the 5.2 ka horizon scale with Mn (Co/Mn = 0.022, Ni/Mn = 0.039;
  Spearman ρ(Co,Mn) = 0.81, ρ(Ni,Mn) = 0.75), while Cu rises ×4 and Ce/La
  from ~1.0 to 1.5–2.3.

## Reproduced by

`companion_analysis/crosslab_matrix_check.py` prints all numbers cited in
Supplementary Methods 14.3, 14.4 and Supplementary Figure 18, and re-renders
Supplementary Figure 18 from this file.
