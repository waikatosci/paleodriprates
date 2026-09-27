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

## 23 Sep 2026: documentation, Fig. 5a marker, calibration series

- WITHDRAWN the same day: the photographic composite of the polished section (`companion_analysis/hs4_composite/`,
  `external/HS4_composite.*`, `external/HS4_axial_slice.*`, Supplementary Figures 21-22 and Fig. 5 panel c).
  The photographs belong to C. Hu and are the subject of his own manuscript on growth rate and mean annual
  temperature; he does not consent to their use here. He also notes that the trace-element powders were not
  drilled on the photographed face (they follow the d18O and d13C profile), so sample depths cannot be located
  on it and the argument built on that placement does not stand. All of it is removed from the repository, the
  manuscript and the Supplementary Information, and nothing else in the paper depended on it.
- README revised: repository tree brought up to date (calibration/, dripwater/, companion_analysis/,
  extended_data/, supplementary_figures/, PRODUCTION_SETTINGS.md), the HS4 inputs explained (full
  589-row table against the 568-sample canonical input), optional dependencies listed, and the AI
  statement updated.
- `dr_app/QUICKSTART.md` now points to `HS4_TE_canonical.csv` (568 samples) rather than the full
  `HS4_TE.csv`, and `dr_app/PRODUCTION_SETTINGS.md` is revised to rev 4: the 568-sample input (the
  retired rev 3 named a 588-row file that kept the second-laboratory samples), acceptance checks and
  point counts updated (568 native; 566 within the dated span).
- Fig. 5a: the modern-baseflow star (17 drips/min at -0.060 ka) sat outside the axis, which ended at
  -0.045 ka, so it was clipped by the right spine and effectively invisible. The young end of the
  axis now extends past the marker (`figures/Fig5_record.py`); nothing else in the figure changes.
- MS text: the one-step precipitation calibration regresses drip rate on the HESHANG P/T series
  (workbook sheet 05b_calibration, PT_h; Supp. Methods 5), not Yichang. Three places in the main
  text said Yichang; corrected. Yichang remains the independent level check in Fig. 4c.

## 24 Sep 2026: growth seasonality (Supplementary Methods 16)

- New `companion_analysis/dripwater_integration_timescale.py`: asks the 2007-2015 monitoring record over
  what window the dripwater Me/Ca ratios track drip rate. The correlation rises monotonically with the
  averaging window (Co/Ca 0.30 at the sampling date, 0.42 at 90 d, 0.64 at 1 yr, 0.73 at 2 yr; Ni/Ca 0.19,
  0.22, 0.36, 0.45), and the seasonal cycle is damped from 0.67 of the mean in drip rate to 0.27 (Ni/Ca),
  0.11 (Co/Ca) and 0.09 (Ca). Summer (Jun-Sep) averages 18.6 drips/min against 13.1 in winter, a factor of
  1.4, so with mu anchored to the annual baseflow mean (16.66) fully summer-weighted calcite carries a
  level offset of ~12%. Output `manuscript_figures/output/TableS_dripwater_integration.csv`.
- New Supplementary Methods 16 uses these numbers to attribute the proxy to baseflow despite the summer
  growth bias: the supply is integrated upstream, the seasonal contrast in flow is small against the
  calibration anchor, event amplitudes are internal to the record and cancel a constant weighting, and any
  residual bias reads wetter, so a drought is understated rather than exaggerated.

## 24 Sep 2026: Fig. 5 legibility and the Discussion restructure (B. Fox review)

- `figures/Fig5_record.py`: the median and interquartile lines over panel a are now white with a
  dark outline (`path_effects.withStroke` in `COL_MEDIAN`) at both native and 1-cm resolution. Over
  the dark end of the posterior colourmap the previous dark lines were lost; the outlined white
  lines read against both ends of the scale.
- The modern-baseflow star is now drawn beside the left-hand drip-rate axis (a blended transform
  with a 7-point outward offset and `clip_on=False`) rather than inside the axes at -0.060 ka. The
  axis extension added on 23 Sep to keep the marker visible is reverted, so the x-axis again ends at
  the youngest sample. The marker is a level annotation, not a datum, and belongs on the axis.
- The 230Th tie-points and their 2-sigma bars were already drawn above panel a; the Fig. 5 caption
  now says so.
- Manuscript only (no code): the Discussion is divided into "Climate reconstruction", "Reconciling
  the drip-rate and magnetic records", "Limitations" and "Implications for hydrological
  reconstruction", with the limitations paragraph moved below the magnetic-record section. The
  passage arguing from a "two-stage regression" is withdrawn: the conversion has used a single
  regression since round 3, and the compression of reconstructed extremes is restated in the
  Limitations section in terms of the intercept of that one regression and of age propagation.
  The drip-rate-to-rainfall arithmetic is stated once, in the Methods, and the Results now point
  there. Reference numbering must be checked at the end, since the moves change citation order.

## 24 Sep 2026: Discussion reference audit (manuscript only, no code change)

Every attribution in the Discussion was paired with the work cited. Corrections applied:

- Mayewski et al. 2004 was credited in the Methods with the rolling z-score used for soft-magnetic
  event detection. That paper screened ~50 records against the Denton and Karlen glacier framework
  and used no such method. The attribution is withdrawn and the reference moved to the 5.2 ka
  paragraph, where the 6.0-5.0 ka interval of globally expressed rapid climate change is what it
  actually supports.
- The claim that speleothem d18O integrates upstream rainout, moisture source and seasonality was
  cited to a vegetation-feedback paper; it now cites Lachniet (2009), which is that paper's subject.
- Aerosol-driven suppression of the EASM was cited to Moon et al. (2023), which attributes
  intensification to greenhouse warming and explicitly not to aerosols, and to Lu et al. (2025), a
  palaeo tipping-point study with no modern large ensemble. Both are removed from that sentence,
  which now cites Lang et al. (2025) alone. Lu et al. moves to the sentence on the reconstructed
  Holocene range, which is what it provides.
- The projection references were split from the reconstruction reference so each sits with the claim
  it supports.
- Episodic delivery of fine magnetic minerals now cites Zhu et al. alongside the karst flow-path
  reference, and the Figure 7 caption carries the numeral for Zhu et al. rather than the name alone.

Reference NUMBERING is not corrected here and should not be: the citations are being relinked to
Zotero, and a refresh regenerates the sequence and the reference list from the field codes.

## 24 Sep 2026: S. Breitenbach's readability edits merged (manuscript only)

His copy branched from the version circulated before the Heshang correction and before B. Fox's
round, so it was merged sentence by sentence rather than taken wholesale.

Taken: the sentence splits, which are the substance of his review. The magnetic-record paragraph
now runs as three sentences rather than one colon-joined chain; the forcing clause in the 5.2 ka
paragraph is its own sentence; the d13C sentence is broken in two; "The synchrony of this
excursion with" replaces "Its synchrony with"; "broader reorganisation of (sub)tropical
atmospheric circulation"; "the dominant source of uncertainty in the quantification of past
precipitation is the error in palaeotemperature estimates"; "stalagmite growth rate" for "calcite
extension"; "prolonged, stable climate states"; "At large" dropped; the model-validation clause
lifted out of its participial tail. The limitations paragraph is split in two, between the
conversion caveats and the calibration assumption, and "This buffering" is renamed "Epikarst
buffering" so it does not depend on an antecedent across the break.

Not taken: his opening to the magnetic section, which is the pre-Fox "Beyond reconciling the two
archives" and reintroduces the dangling reference she flagged; his restored forward reference to
the magnetic record in the limitations; his citations at 5.2 ka, which predate this morning's
corrections; his Methods paragraph, which predates the Heshang correction; and "5.2 ka BP"
throughout, since the figures and Supplement use "ka".

CORRECTED: he had read the "~70%" at 5.2 ka as the fraction of percolation remaining rather than
the size of the fall, and edited the sentence to say percolation "remained above two-thirds" of
its mid-Holocene level, describing a "ca. 30% reduction". That inverts the central result. The
sentence now reads "dropped to roughly a third of its mid-Holocene level, a reduction of about
70%", which states both and removes the ambiguity he read past.

work/validate_reject2.py added: the earlier check counted a tracked paragraph split as a new
paragraph and failed a correct document. The new one merges paragraphs whose mark is an insertion,
which is what rejecting the split actually does.

## 24 Sep 2026 (later): findings resolved against the data; S. Breitenbach's remaining edits and comments

Resolved against the reconstruction outputs (Fig6_P_recon_dense, Fig6_P_comparison.py, the
native-resolution summary with 02_chronology ages):

- Fig. 6b caption said the 5.2 ka window shows a "~28% negative anomaly". The script computes each
  bar as the window median against the 0-9.5 ka Holocene mean of the dense series (1,107 mm/yr);
  the 5,200-5,100 window median is 897 mm/yr, -19%, and that is what the rendered panel shows.
  The 28% (26-27% on the data) is the departure from the 6.5-5.5 ka mean. Caption corrected to
  ~19%; the Methods sentence on baselines now states what the script does.
- The 5.2 ka precipitation minimum (~896 mm/yr) equals the late-Holocene (2.0-0.2 ka) mean
  (905-918 mm/yr). That is a real feature of the reconstruction, not an error, and the manuscript
  should say so rather than let a referee find it; wording pending the first author.
- "The post-1800 surge" in the Limitations is the 1770-1850 CE surge named everywhere else: onset
  between the 1767 (6.3 drips/min) and 1777 (41) samples, peak 84 at 1805 CE. Label aligned.
- The post-1850 decline (-0.04 drips/min/yr, R^2 = 0.70) is the native-series regression from 1850
  CE to the top of the record (1995 CE; n = 34). The text now gives the window. A referee reading
  the slope against the 1770-1850 median of 26 would get 225 years; the fit starts at 22 in 1850
  and reaches 17.6 by 1960, which is what the sentence says.

S. Breitenbach's edits outside the four Discussion paragraphs merged sentence by sentence (44
tracked changes, two paragraph splits: the detrital-screening sentences off the teleconnection
paragraph in the Results, and the sample-to-sample-variability sentences off the recent-record
paragraph in the Discussion). Not taken: "ka BP" throughout; "novel"; "median chronology" for
"best-estimate chronology"; "conservative upper bound" (it is the opposite); "triggered" for
"contributed to" on the Nile Valley; his split of the delta-18O sentence, which made the
integration explain the sign rather than the weak coupling; "strong" for "coherent" floods;
"Salmanzadeh" without the initial. His 38 comments are carried across with his authorship and
dates, anchored to the passage each concerns or its rewritten successor.

Toolkit: the helper regex for tracked deletions matched self-closing paragraph-mark deletions
(<w:del .../>) and swallowed the live text up to the next </w:del>; fixed in the r19 script
(`<w:del [^>]*[^/]>`), and the same fix belongs in edits.py if the helpers are reused.

## 24 Sep 2026 (later): numbers reconciled against the workbook; remaining S. Breitenbach edits

Manuscript only. Four statements the independent read flagged were checked against
`HS4_SourceData.xlsx` (Fig6_P_recon_dense) and the drip-rate summaries on the 02_chronology age model:

- 5.2 ka precipitation anomaly. The 5,200-5,100 window median is 898 mm/yr against a 0-9.5 ka mean of
  1,107: -19%, which is what Fig. 6b draws. The Fig. 6b caption said ~28% (stale) and the Methods
  described baselines the figure does not use (late-Holocene mean; mid-Holocene mean for the event).
  Caption now ~19%; Methods now state the Holocene-mean baseline and give the event relative to the
  surrounding mid-Holocene mean (6.5-5.5 ka, 1,216 mm/yr: -26%, "roughly a quarter").
- ~900 mm/yr is also the late-Holocene (2.0-0.2 ka) mean (918). Not an error: the event drops to a
  level the record does not otherwise reach until the late Holocene. Fig. 6a caption now says so.
- Surge dates. Native series: first sample above 1.5x the late-Holocene median at 1777 CE (text
  "~1770 CE", kept); peak 83.7 drips/min at 1805 CE ("the post-1800 surge" now "the surge, which peaks
  ... around 1805 CE"); the figure band "post-1750 CE" is a round-number band and is kept.
- Post-1850 decline. Native median regression 1852-1995: -0.044/yr, R2 0.70, as stated; it starts at
  ~21 drips/min (1850-1875 median 20.9), not at the 1770-1850 median of 26, and the first sample at or
  below 17.5 is 1952 CE. The text now gives the starting level, so the arithmetic closes.
- The "+10%" post-1800 precipitation anomaly is the 1750-1850 window median against the late-Holocene
  mean (+9%); the baseline is now named.

S. Breitenbach's edits outside the four Discussion paragraphs, taken where correct and better: the
Abstract ("site-scale hydrological dynamics"; "Holocene" dropped from the closing sentence), the DGT
sentence split, "highly consistent", "at this drip point", the Fig. 4c caption, the orbital-decline
paragraph (units on the age-propagated range, "This decline mirrors" as its own sentence, "significant"
on the p < 0.01 correlation), "background of changing orbital forcing", the conservative-figure
sentence at 5.2 ka, "compared with", the second-laboratory sentence split, the detrital exclusion as its
own paragraph, "In marked contrast to the 5.2 ka event", "event-specific precipitation anomaly", "not
significantly reduced", the age-inversion sentence split, "PCP signature", past tense at 8.2 ka, the
Qing-drought and 1770-1850 sentences split, and the resolution remark as its own paragraph. Not
taken: "5.2 ka BP" throughout, his Methods paragraph (predates the Heshang correction), his edit to
affiliation 9 (C. Hu's, left as supplied), and his unfinished "xxx" sentence.

His 38 comments moved into the live file (37; one was empty), re-anchored on the current text: 24
exact, 7 disambiguated by context, 5 on the surviving fragment, 1 by preceding context; nine placed
by hand where his anchor text no longer exists. Comment ids offset by 2000; his paraIds and durableIds
kept so the Extended/Ids/Extensible parts carry across. Accepted text verified byte-identical before
and after the transplant.

## 24 Sep 2026: Cr at 5.2 ka (manuscript only)

The Results gave Cr at the event as "at most ×1.8" in one paragraph and "×1.36, p = 0.06" in the
next. Both reproduce: ×1.8 is EF_Cr_max in the 5.2 ka core against the local rolling baseline
(detrital_ternary_screen.py; TableS_event_enrichment), ×1.36 is the window-median contrast against
the rest of the record (events_differ_lithogenic.py, S13). The text now says which is which.

## 24 Sep 2026: TraCE product name (manuscript only)

The Fig. 6 comparator is TraCE-21K-II (sheet Fig6_TraCE, from
East_Asia_TraCE21KII_PRECT_mm_yr_LGM_to_present.csv). The Methods called it TraCE-21ka and both
mentions cited only Liu et al. 2009. Both now read TraCE-21K-II and cite Liu et al. 2009 with He &
Clark 2022 (Nat. Clim. Change 12, 449-454), added as reference 58.

## 24 Sep 2026: PCP expansion (manuscript only)

The Abstract expanded PCP as "prior carbonate precipitation" and the Introduction as "prior calcite
precipitation". The Abstract now matches the Introduction: the system precipitates calcite, and
"carbonate" would leave aragonite open.

## 24 Sep 2026: Fig. 5 legibility changes reverted

The two `figures/Fig5_record.py` changes made after the B. Fox review (white outlined median and
interquartile lines; modern-baseflow star drawn beside the axis) are withdrawn at the lead author's
request: the white lines are not an improvement, and the star had been placed at the old end of the
axis. The script returns to the 23 Sep version: dark median and interquartile lines, and the star at
its own age (-0.060 ka, 17 drips/min) at the young end of panel a, with the x-axis extended past it
so it sits inside the frame. Outputs re-rendered.

## 24 Sep 2026: Discussion streamlined (manuscript only)

On the lead author's returned working copy: the Discussion tightened for readability at the
co-reviewers' request (repetition of Results figures removed, long sentences split, the magnetic
section reduced to two paragraphs), the Fig. 5 caption and image brought back to the reverted
render, the Fig. 6 caption colours matched to the figure (HS4 blue, TraCE-21K-II red), and the
8.2 ka flanking-test citation pointed at Supplementary Fig. 14. No code or data changed.

## 25 Sep 2026: editorial pass (manuscript only)

Results and Discussion wording tightened for precision (kinetic, not thermodynamic; metal–organic,
not organometallic; interpretation moved out of the Results), and growth continuity through the
5.2 ka and 8.2 ka intervals stated from the chronology (02_chronology tie-points at 154.0 and
165.8 cm give 0.32 mm/yr across the 5.2 ka bracket; annual laminae are counted through 8.2 ka by
Liu et al. 2013). No code or data changed.

## 28 Sep 2026: release metadata for v1.1.0

Added `CITATION.cff`, `.zenodo.json` (title, three code authors, MIT, QUEST grant, links to the
GitHub repository and to the data deposit 10.5281/zenodo.16392750), `RELEASE_NOTES_v1.1.0.md` and
`DATA_DEPOSIT_MANIFEST.md` (the files for the new version of the data record, and the two that the
repository does not hold). No code or data changed.
