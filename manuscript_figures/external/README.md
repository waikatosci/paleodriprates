# external/ — inputs not stored in the workbook

Some inputs are too large for a spreadsheet or are third-party published products. They live
here, each with its role, size and sha256. Where a value here feeds a plotted number, the
derivation is recorded on the relevant workbook sheet.

## In repo

| file | used by | note |
|---|---|---|
| `pdf_heatmap_hr.json` | Figs 5a, 7 | drip-rate PDF heatmap (210 V-bins x 588 depths) from the canonical native run re-executed headlessly on 2026-09-16 with the PRODUCTION_SETTINGS parameters (input summary in `run_logs/input_summary_native.csv`); summary percentiles identical to `drip_rate_summary_hr_censored.csv` |
| `pdf_heatmap_lr.json` | (reference) | 1 cm run heatmap, edge nodes trimmed to the 255 data-spanning nodes; percentiles identical to workbook 04_driprate_posterior (`run_logs/input_summary_1cm.csv`) |
| `drip_rate_realisations_hr.csv.gz` | RQA (`companion_analysis/RQA_HS4_ensemble.py`) | 1,000 MC realisations of the native run (seed 42) |
| `drip_rate_summary_hr.csv` | Figs 5a, 7 | native-resolution drip percentiles, canonical Run 1 (2026-07-20; mu canonical, sigma = pi/sqrt(6), lambda_F = 0.01): 586 pts, the two below-resolution points (156.23, 156.64 cm) dropped; full 588-row table with censoring flags in `drip_rate_summary_hr_censored.csv` |
| `drip_rate_summary_lr.csv` | Fig 7 | 1 cm-interpolated drip percentiles, canonical Run 2 (255 pts, edge nodes trimmed); identical to workbook 04_driprate_posterior |
| `HS4_Zhu2017_IRMsoft_flux.csv` | Fig 7 | IRM_soft flux (Zhu 2017), 115 samples; columns age_mid_kaBP, IRMsoft_flux_Am2_per_yr |

## Fetch at runtime (third-party)

| file | used by | source |
|---|---|---|
| `precip.mon.mean.nc` | Fig 2a/2b | GPCP v2.3 monthly precipitation (NOAA PSL) |
| `SR_LR.tif` | Fig 2a | Natural Earth 10m shaded relief |

Retrieval URLs and checksums are in `fetch_external.sh`.

## Added for the revision (companion-analysis inputs)

| file | used by | note |
|---|---|---|
| `HS4_TE_multielement.csv` | `detrital_ternary_screen.py`, `make_fig_events_differ.py` | full solid-phase suite (V, Cr, Co, Ni, Cu, Zn; 589 rows) exported from workbook `01_master_TE` by `calibration/export_multielement_te.py` (single documented 0.06 cm exclusion applied) |

The two source-stability companion scripts read `dripwater/HS4_dripwater_canonical.csv`
(92 paired dripwater analyses, 2007-2016) at the repo root, not this folder.

| `HS4_8p2ka_digest_2017.xlsx` | SM13.2–13.3 (`companion_analysis/events_differ_lithogenic.py`) | full-suite carbonate-digest ICP-MS, 51 samples across 230.8–236.8 cm (30 May 2017); also copied into workbook sheet 07_lithogenic_8p2ka |

2026-09-17 RE-ANCHORING: all drip-rate products above regenerated at mu_Ni = -3.8372, mu_Co = -5.4171 (calibration target 16.66 drips/min = 2004-2023 annual-baseflow mean at V_DROP 0.14 mL; run logs in run_logs/). Level +10.4% vs the 15.09-target run, shape unchanged. `drip_rate_summary_ap.csv` = age-propagated run (1,000 realisations, 1-yr grid, V grid 2,000; = workbook Fig5_driprate_AP, not rescaled). `calibration_onestep.csv` and `T_recon_Wang_et_al.xlsx` feed manuscript_figures/build_precip_onestep.py (Fig 6 sheets).
| `drip_rate_realisations_ap.csv.gz` | `drip_rate_stationarity_tests.py` (SM7), `companion_analysis/RQA_HS4_ensemble.py` (SM6, Supp Fig 12) | 1,000 MC realisations of the age-propagated run (seed 42; 10,293 annual steps; ~100 MB, NOT tracked in git — regenerate with Dr Paleo in age mode at the PRODUCTION_SETTINGS parameters, or download from the Zenodo deposit); full stationarity-test output in `run_logs/stationarity_tests_ap.log` |

2026-09-21 SECOND-LABORATORY EXCLUSION: all drip-rate products above regenerated on the 568-point canonical input (`dr_app/HS4_example_inputs/HS4_TE_canonical.csv`) after removal of the 20 samples from the supplementary analytical run at a second laboratory (12 HS4-A- at 155.2-157.8 cm; 8 HS4-C- at 7.0-8.5 cm; listed with reasons in `calibration/excluded_points.csv`). The data owner (C. Hu) reports that the second-laboratory calibration does not reproduce a certified standard (Ni 5.20 ppm returned 2.38 ppm); at the three depths analysed by both laboratories (8.09, 156.38, 157.48 cm) the second-laboratory values are linearly related to the primary values (r >= 0.998) with an offset of ~ +3.9 ppm Ni and ~ +0.9 ppm Co. No point in the cleaned record falls below the joint resolution limit, so `drip_rate_summary_hr.csv` and `drip_rate_summary_hr_censored.csv` now carry the same 568 rows (censored = 0 throughout). Native run: 24 s; 1 cm run: 19 s; age-propagated run: 618 s (1,000 realisations, seed 42). Run logs in run_logs/. The 588-point products of 2026-09-17 are retained outside the repository for comparison.
