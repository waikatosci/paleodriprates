# Round-2 revision patch — second-laboratory full ICP-MS suite

Add-on to the paleodriprates repository at commit `fb3fc80`.

## Files to add (new)

| path | size | sha256 |
|---|---|---|
| `manuscript_figures/external/HS4_TE_full_suite_both_labs.xlsx` | 51,375 B | `72dd17205b3fd0b133186e41dad8cb213a4c7d18003ff2fb57f79b061750b1e7` |
| `manuscript_figures/external/HS4_TE_full_suite_both_labs.README.md` | 2,476 B | (this note) |
| `companion_analysis/crosslab_matrix_check.py` | 12,508 B | (this note) |

## Files to update (edit in place)

- `calibration/excluded_points.csv` — reason field on the 20 second-lab rows
  updated from "certified-standard failure" to the Ca-matrix explanation
  (the corrected file is `excluded_points.csv` in the delivery bundle).
- `manuscript_figures/external/README.md` — add one row under the
  "Added for the revision" table:

  ```
  | `HS4_TE_full_suite_both_labs.xlsx` | `companion_analysis/crosslab_matrix_check.py`, Supp. Methods 14.3–14.4, Supp. Fig. 18 | full trace-element suite (Li–U) for 20 HS4 powders (plus a procedural blank) and four silicate reference materials; the 2019 second-laboratory run that supplied the twenty samples now withdrawn from the input, plus the paired 2009 primary-run values at the three shared depths (see companion README) |
  ```

- `companion_analysis/crosslab_sensitivity.py` — replace the docstring
  reference to a "certified-standard failure" with a pointer to
  `crosslab_matrix_check.py` and the Ca-matrix interference. The
  sensitivity analysis and its numbers do not change.

## Commit message

```
Add the second-laboratory full ICP-MS suite and identify the Ca matrix as the cause of its offset

The second laboratory (Wuhan, November 2019) measured the full trace-element
suite for 20 HS4 powders (12 across the 5.2 ka event, 155.2-157.8 cm; 8 near
the growth surface, 7.0-8.5 cm) and returned AGV-2, BHVO-2 and BCR-2 within a
few percent of certified values. It nonetheless offsets from the primary run
at the three depths analysed by both laboratories. The cause is the calcium
matrix: outside the 5.2 ka core, Co, Ni and Fe in the second run are linear
functions of Ca in the digest (r = 0.98, 0.96, 0.97) with a contribution at
calcite Ca of ~1.6 ppm Co, ~5.8 ppm Ni and ~6,900 ppm Fe (~0.7%). The Fe is
the 40Ca16O+ polyatomic interference at m/z 56; the Co and Ni baselines are
consistent with the corresponding Ca-bearing ions at m/z 59 and 60. The
primary run's standards were matrix-matched to high-purity Ca and do not
carry it. The 20 samples remain withdrawn from the inversion for that
reason.

The same file also provides Al, Th, Ti, Zr and Mn at the 5.2 ka horizon,
which the primary run did not measure. Al up to 38 ppm, Th up to 0.013 ppm
and Ti up to 7 ppm bound any detrital contribution at a few per cent of the
observed Co and Ni excesses at most (1-2 percent for the largest excesses); the excess Co and Ni scale with Mn (Co/Mn 0.022,
Ni/Mn 0.039) with no lithogenic response, the signature of a residence-time
excursion in four organically bound metals. This is the mass balance quoted
in Supplementary Methods 14.4 and Supplementary Figure 18.

  add:   manuscript_figures/external/HS4_TE_full_suite_both_labs.xlsx
  add:   manuscript_figures/external/HS4_TE_full_suite_both_labs.README.md
  add:   companion_analysis/crosslab_matrix_check.py
  edit:  manuscript_figures/external/README.md
  edit:  calibration/excluded_points.csv
  edit:  companion_analysis/crosslab_sensitivity.py    # docstring only

Reproduces every value in Supplementary Methods 14.3, 14.4 and Supplementary
Figure 18 by running `python companion_analysis/crosslab_matrix_check.py`.

Data provided by C. Hu, 21 September 2026.
```

---

# Round-3 addendum (22 September 2026)

Applied on top of the round-2 patch; every number below is printed by
`companion_analysis/crosslab_matrix_check.py` from
`manuscript_figures/external/HS4_TE_full_suite_both_labs.xlsx`.

## Corrections

- The 2019 run holds **20** HS4 powders, not 21: the 21st column of sheet
  `Lab 2` is `HS1-0`, a procedural blank. Count corrected in the script,
  both READMEs, this note and Supplementary Methods 14.4 / Figure 18.
- The detrital share of the Co and Ni excess is **1.0–3.5% (Co) and
  1.6–6.0% (Ni) sample by sample**, 1–2% for the samples with the largest
  excesses; "below 1%" was not supported by the numbers as printed.
- Spearman ρ between the excess Co and Ni (over the Ca-proportional
  baseline) and Mn is **0.85 and 0.70** (0.81 and 0.75 were the raw-
  concentration values).
- The light rare earths are enriched roughly tenfold in the 5.2 ka core,
  but Ce/La in the core (1.0–2.3) overlaps its range outside it (0.8–1.7)
  and does not track Mn: **no cerium anomaly is resolved**. The two samples
  formerly called "the two Mn-richest" (156.23 and 156.64 cm) are the two
  with the highest Co and Ni; by Mn they rank first and third.
- A-880 (157.01 cm, the resolved 5.2 ka minimum) has Th 0.010 ppm and
  Zr 0.35 ppm in the **primary run's own full-suite analysis** (sheet
  `lab 1`, W0907041), not in the 2019 run.
- `detrital_ternary_screen.py` now takes V_MODERN / V_FLOOR (15.62 / 2.81)
  from `source_variation_propagation.py` instead of the pre-refit 14.14 /
  1.07; the kinetic-path slope quoted in Supplementary Figure 8 is unchanged
  (1.22, "~1.2").
- `crosslab_sensitivity.py` docstring: the certified-standard explanation
  replaced by the Ca-matrix one (the sensitivity result is unchanged); the
  docstring patch file is removed.

## Additions

- Leverage test: without the highest-Ca sample (HS4-A-873, 536 ×10³ ppm)
  the Ca correlations remain significant (r = 0.88 Co, 0.68 Ni, 0.80 Fe;
  p ≤ 0.014, n = 12); Co and Fe slopes unchanged, Ni slope 4.2 ppm at
  calcite Ca. Reported Ca runs 376–536 ×10³ ppm (nine samples above the
  stoichiometric 400 ×10³ ppm), i.e. an analytical quantity of each digest.
- Supplementary Figure 18 re-rendered with the panel-a legend above the
  data (the previous position hid two of the seven core samples).
- `extended_data/README.md` with the working invocation of each SF11–16
  script; `astropy` added to `requirements.txt` for SF15.

## Merged from the July revision bundle (nc-monsoon-code, 22 Sep 2026)

`dripwater/` (SF6, SF7, SM11.6), `companion_analysis/events_differ_lithogenic.py`
and `make_fig_events_differ.py` (SM13, SF9; the digest now reads from workbook
sheet `07_lithogenic_8p2ka` when `HS4_8p2ka_digest_2017.xlsx` is absent),
`age_model.csv` in each `extended_data/` folder (SF11, SF13, SF16),
`extended_data/Ex_Data_2_RQA/` (SF12), `dr_app/PRODUCTION_SETTINGS.md`,
`COMMIT_PLAN.md`. `manuscript_figures/external/calibration_onestep.csv` was
rebuilt from workbook sheet `05b_calibration` (identical values) and
`T_recon_Wang_et_al.xlsx` copied to `external/`, so `build_precip_onestep.py`
runs (peak 1,402 mm, late 922 mm, decline 34%, R² 0.303, as in the text).
`Ex_D1.py` defaults now carry the 2004–2023 monitoring values of the SF11
caption (CV 0.445 raw, 0.027 corrected); `ExD_3.py` finds `Drip_rate.xlsx` at
the repository root; `make_fig_events_differ.py` writes its 2σ annotation from
the data (200 vs 60 yr) instead of a stale hard-coded 210 vs 72.

Reproduced and checked against the text on 22 Sep 2026: SM11.6/SF7
(±35% mid-record, ±39% at the floor, ≤20% of the range), SM12/SF8, SM13/SF9
(ρ(Al,Co) +0.31, Mg/Ca–Sr/Ca +0.84, Cr ×1.18 p 0.002, δ¹³C +1.63‰, coupling
+0.05 / +0.89), SM14.3–14.4/SF18, SF11, SF13, SF14, SF15, SF16, Fig 6 series.
One correction to the text came out of it: the V enrichment at 5.2 ka is
×1.8 with p = 0.22 on the four primary-run samples (S13.3 said ×2.0,
p = 0.001, a value from the 589-row table that still held the withdrawn
samples).

## Supplementary Figure 17 made reproducible

`companion_analysis/run_crosslab_585.py` drives Dr Paleo headlessly at the
logged production settings. `--check` re-runs the 568-point input and
reproduces `drip_rate_summary_hr.csv` exactly (568/568 depths, 0.000%);
the default builds the 585-point input (568 primary + 17 rescaled
second-laboratory rows; fits Ni 1.0838x − 3.9022, Co 0.9195x − 0.8892) and
writes `external/drip_rate_summary_hr_corr585.csv` plus its parameter log.
Minimum 1.61 drips/min at 156.23 cm, 91% below the surrounding median of
18.6, against 2.81 (85%) for the record used: the values in S14.3 and the
SF17 caption. `crosslab_sensitivity.py` now reads the released summaries
(no private run folders).

## Files the repository still needs

- `manuscript_figures/external/drip_rate_realisations_ap.csv.gz` (SM7, SF12;
  regenerate with Dr Paleo in age mode, or take from the Zenodo deposit)
- the DGT deployment table behind Supplementary Figure 19, for Source Data

## Round 4 (22 Sep 2026): Ni-constrained 5.2 ka minimum

- New `companion_analysis/run_single_metal.py`: the canonical inversion on Ni alone or Co alone
  (native, 1-cm and age-propagated modes; `--metal both` re-runs the joint inversion as a check
  and reproduces the released summaries to 0.000%). Outputs
  `manuscript_figures/external/drip_rate_summary_{hr,lr,ap}_Nionly.csv`,
  `drip_rate_summary_hr_Coonly.csv` and their run logs.
- New `companion_analysis/ni_co_agreement.py`: Supplementary Figure 19 and the numbers of
  Supplementary Methods 14.4. Outside the 5.2 ka interval the Ni- and Co-implied drip rates agree
  (ratio median 1.15, 5th-95th percentile 0.80-1.51, n = 565); at the three 5.2 ka samples Co gives
  drip rates about three times lower than Ni, the three largest departures in the record. The
  event is now reported from Ni alone: minimum 6.42 drips/min at 157.01 cm (IQR 6.11-6.75), 69%
  below the surrounding median (71% against the 5,500-6,500 BP background; 6.46 and 71% at 1 cm;
  11.1 drips/min and ~50% age-propagated). The joint minimum (2.81, 85%) is kept as an upper bound.
- `run_crosslab_585.py`: `production_params()` takes the run log to read, and `run()` accepts
  extra input files (the age-depth table for the age-propagated mode). Behaviour of the script
  itself is unchanged.
- The DGT figure is now Supplementary Figure 20.

## Round 5 (22 Sep 2026, afternoon): second-laboratory run as replication

- New `companion_analysis/crosslab_replication.py`: corrects each 2019 value for the Ca of its own
  digest (slope from the 13 non-core samples, relative to the run mean), places it on the
  primary-run scale with the relation at the three shared depths (Ni: primary = 1.117 x corrected
  second - 4.222), and runs the Ni-only inversion on the 568 + 17 samples. The corrected second-run
  samples reproduce the 5.2 ka event on Ni (5.1-5.9 drips/min at 156.23-156.88 cm beside 6.4 at
  157.01 cm; nine samples in 155.2-157.0 cm; minimum 5.12, 75% below the surrounding median). The
  Ca step removes the false dip at 157.65 cm (the Ca-richest digest; 9.4 with the paired-depth
  relation alone, 18.4 with the Ca step). Outputs `drip_rate_summary_hr_Nionly_rep585.csv`,
  `crosslab_replication_2019_corrected.csv` and the run log. The primary record is unchanged.
- `ni_co_agreement.py`: Supplementary Figure 19b now overlays the corrected second-run Ni values.
- New `companion_analysis/element_covariance.py`: which elements covary with Ni and Co in the
  record, the 2019 full suite and the dripwater (tables and a summary figure in
  `manuscript_figures/output/`; not cited in the paper).
- Methods: the primary (2009) run did not measure Ca (determined separately by ICP-AES, per
  C. Hu, 22 Sep 2026); the earlier statement that its standards were matrix-matched with Ca is
  withdrawn.

## 23 Sep 2026: second-laboratory details (C. Hu)

- The 2019 run was made at Wuhan SampleSolution Analytical Technology Co., Ltd on an Agilent 7700e
  ICP-MS with the method of the 2009 run (not an Agilent 8900 triple-quadrupole, as previously
  stated). Docstring of crosslab_matrix_check.py and the full-suite README corrected.
- C. Hu has agreed to the release of the measurement report as Source Data.

## 23 Sep 2026: photographic composite of the section

- New `companion_analysis/hs4_composite/`: `build_composite.py` stitches C. Hu's eight photographs of the
  polished section onto the pencil depth scale (100 px/cm; tick residuals < 1 mm; overlaps agree to
  0.1–0.5 mm) and cuts a 6 cm axial slice; `si_figures.py` draws Supplementary Figures 21 (composite) and
  22 (the 5.2 ka horizon with Ni and Co). New inputs `external/HS4_composite.*` and `external/HS4_axial_slice.*`.
- The enriched 5.2 ka samples lie within the loose fragment between the breaks at 155.1 and 157.3 cm, with the
  highest Ni and Co in its interior (156.2–157.0 cm) and background values at its lower break (157.10–157.30 cm).
- Fig. 5 gains panel c: the axial slice placed on the age axis through the best-estimate age model, depth marked
  below (`figures/Fig5_record.py`; the event bands carry through).
- Band spacing in the photographs was tested against the 230Th growth rate and does not track it; the
  photographs do not resolve annual layers, so no lamina counts are reported (README in the new folder).
- README revised: repository tree brought up to date (calibration/, dripwater/, companion_analysis/,
  extended_data/, supplementary_figures/, PRODUCTION_SETTINGS.md), the HS4 inputs explained (full
  589-row table against the 568-sample canonical input), optional dependencies listed, NSFC funding
  and the photographs acknowledged, and the AI statement updated. `opencv-python-headless` added to
  requirements.txt for the composite.
- `dr_app/QUICKSTART.md` now points to `HS4_TE_canonical.csv` (568 samples) rather than the full
  `HS4_TE.csv`, and `dr_app/PRODUCTION_SETTINGS.md` is revised to rev 4: the 568-sample input (the
  retired rev 3 named a 588-row file that kept the second-laboratory samples), acceptance checks and
  point counts updated (568 native; 566 within the dated span).
