# external/ — inputs not stored in the workbook

Some inputs are too large for a spreadsheet or are third-party published products. They live
here, each with its role, size and sha256. Where a value here feeds a plotted number, the
derivation is recorded on the relevant workbook sheet.

## In repo

| file | used by | note |
|---|---|---|
| `pdf_heatmap_hr.json` | Figs 5a, 7 | drip-rate PDF heatmap, 209x576, sigma = pi/sqrt(6) |
| `drip_rate_summary_hr.csv` | Figs 5a, 7 | native-resolution (576 pt) drip percentiles |
| `drip_rate_summary_lr.csv` | Fig 7 | 1 cm-interpolated (258 pt) drip percentiles; Fig 5 reads its copy from workbook 04_driprate_posterior |
| `HS4_Zhu2017_IRMsoft_flux.csv` | Fig 7 | IRM_soft flux (Zhu 2017), 115 samples; columns age_mid_kaBP, IRMsoft_flux_Am2_per_yr |

## Fetch at runtime (third-party)

| file | used by | source |
|---|---|---|
| `precip.mon.mean.nc` | Fig 2a/2b | GPCP v2.3 monthly precipitation (NOAA PSL) |
| `SR_LR.tif` | Fig 2a | Natural Earth 10m shaded relief |

Retrieval URLs and checksums are in `fetch_external.sh`.

