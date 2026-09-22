# Dr Paleo: HS4 production run settings (canonical, rev 4, 2026-09-23)

Reproduces the manuscript record at the locked calibration. Any deviation
from this sheet is not the production record.

## Input file (rev 4, 2026-09-22)
Upload **HS4_TE_canonical.csv** (in `HS4_example_inputs/`): the 568 samples of
the primary (2009) ICP-MS run. It is the full raw TE table (`HS4_TE.csv`, 589
rows, value-identical to Drip_rate.xlsx 2.Trace_Elems) minus 21 documented
exclusions: the 0.06 cm surface-cap point (Ni/Co both >2× neighbours; would
invert to a spurious slow-drip point at the top of the record) and the 20
samples of the second-laboratory run (2019), withdrawn on 2026-09-22 and
examined separately (Supp. Methods 14.3–14.4; companion_analysis/crosslab_*.py).
Do not upload `HS4_TE.csv` for a production run; it still holds those rows.

Rev 3 (2026-07-20) used a 588-row input that kept the 20 second-laboratory
samples (file `HS4_TE_canonical_588.csv`, no longer in the repository).
Build chain lives in the repo: calibration/excluded_points.csv (the rule) +
calibration/make_canonical_te_input.py (applies it). To change cleaning,
edit the CSV and re-run the script — never hand-edit data files.

RETIRED: HS4_TE_production_576.csv (rev 2) and the legacy 13-row exclusion
behind the published 576-point external grid. Twelve of those rows had no
recoverable rationale and are re-included; note they all sit above their
local Ni baseline (1.24–2.59×) — if a detrital screen is ever re-instated,
do it via excluded_points.csv. The run's own Hampel/Pareto detector
(window 11, ZERO_TOL 1e-3, PDF_TOL 3) median-replaces spikes; on this input
it flags Ni @ 111.10 and Co @ 8.48, 36.03 cm.

## Preprocessing panel ("Block-average / Sigma-clip / Window size")
- **Do not press Apply & Save.** It overwrites the uploaded TE csv in place;
  after any accidental save, re-upload fresh.
- Block-average: unused. Sigma-clip: unused (outlier handling happens inside
  the run). Window field defaults to 11 in the patched app and doubles as the
  run's detector window — leave it.

## App defaults (2026-07-20 patched app.py)
Kd mode defaults to Literature/manual (canonical ELEM_DEFAULTS survive load);
window 11, v_max 100, v_res 5000, n_realisations 1000, rng_seed 42,
temp 12 °C, global drip 16.66, station Heshang. On an UNPATCHED app: switch
both TE cards to Literature mode before loading data, set window 11, and
verify v_max 100.

## Canonical parameter set
| Parameter | Ni | Co |
|---|---|---|
| Kd_mn (ln k_d) | **−3.8372** | **−5.4171** |
| Kd_sd | 1.282549830161864 (= π/√6, fixed a priori) | same |
| F (fast, λ_F) | 0.01 | 0.01 |
| InertF (λ_I) | 0.10 | 0.40 |
| labile (auto) | 0.89 | 0.59 |
| Kp | 1.1 | 4.4 |
| K_e | 1 | 1 |
| aq conc (ppb) | 4.370005432 | 0.460490861 |

Global: Ca = 67.43731092 ppm · analysis_mode = full · hiatus zones = none.
Irrelevant in full + literature mode: cal_pct, global drip rate, temperature
(only active when Kp = −1 theoretical).

## Abort-early check
input_summary.csv is written at run start. Verify: te1_Kd_mn = −3.8372,
te2_Kd_mn = −5.4171, outlier_win_size = 11, v_max = 100. If te1_Kd_mn
reads −3.844, the auto-calibration clobbered the boxes — kill the run.

## Run matrix
| Run | depth_res | with_age | feeds |
|---|---|---|---|
| 1. Native | native | no | drip_rate_summary_hr.csv, pdf_heatmap_hr.json (Figs 5a, 7), realisations → RQA (ED Figs 1–2) |
| 2. 1 cm | 1 cm | no | drip_rate_summary_lr.csv → workbook 04_driprate_posterior |
| 3. Age-propagated ("Run4") | native | yes | Fig5_driprate_AP (then re-anchor to modern baseflow as before) |

## Acceptance checks (native run)
- **568 points, depth 0.18–255.5 cm**, no duplicate depths (rev 3 had 588,
  with duplicates at 8.09/156.38/157.48 from the second-laboratory run);
  `drip_rate_summary_hr.csv` is this run
- Top-of-record pc50 ≈ 15.6 drips/min at 0.18 cm (15.74 in rev 2)
- On the 576 shared depths: single-digit-percent shift vs the old σ=π/√6
  externals (expected e^Δμ: −2.8% via Ni, +5.7% via Co), NOT a factor of ~2
- The 12 re-included points will read as slow-drip-side excursions (all are
  high-Ni) — expected, not a failure; watch 132.38/140.8 cm near the event
  window when re-deriving event magnitudes
- 5.2 ka minimum pc50 2.81 drips/min at 157.01 cm (joint Ni–Co; the Ni-only
  run, companion_analysis/run_single_metal.py, gives 6.42), not 0

## MS reconciliation (assembly checklist)
- Restate from this pipeline: **568** on the native depth grid; **566** within
  the dated span (age model 0.0–253.0 cm excludes rows at 253.5 and 255.5 from
  age-mode products). The 585-point input is the sensitivity run that adds 17
  Ca-corrected second-laboratory samples to the 568
  (companion_analysis/run_crosslab_585.py), not the production record.
