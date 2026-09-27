# Data deposit manifest for the new version of 10.5281/zenodo.16392750

Files to upload (present in the repository snapshot; sizes and sha256 prefixes as of 28 Sep 2026):

| file | size | sha256 (first 16) |
|---|---|---|
| `manuscript_figures/HS4_SourceData.xlsx` | 1,072,433 B | `1b72c213cee3efce` |
| `manuscript_figures/external/HS4_TE_full_suite_both_labs.xlsx` | 51,375 B | `72dd17205b3fd0b1` |
| `manuscript_figures/external/HS4_TE_full_suite_both_labs.README.md` | 3,874 B | `a3e2c41b46314e5f` |
| `Drip_rate.xlsx` | 7,093,298 B | `37fb6c7a4327c62a` |
| `manuscript_figures/external/drip_rate_realisations_hr.csv.gz` | 5,739,229 B | `101ddfc9762f9b10` |
| `manuscript_figures/external/drip_rate_summary_hr.csv` | 77,223 B | `ddf01ef943e63b89` |
| `manuscript_figures/external/drip_rate_summary_lr.csv` | 34,527 B | `2c0ddc9006778360` |
| `manuscript_figures/external/drip_rate_summary_ap.csv` | 1,394,954 B | `70d0d0789ebd61cc` |
| `manuscript_figures/external/drip_rate_summary_hr_Nionly.csv` | 77,012 B | `4f256972457a3f57` |
| `manuscript_figures/external/crosslab_replication_2019_corrected.csv` | 2,212 B | `07f15e875f63ed01` |
| `manuscript_figures/external/HS4_TE_multielement.csv` | 39,289 B | `03e8eedcb55184e2` |
| `manuscript_figures/external/HS4_Zhu2017_IRMsoft_flux.csv` | 4,359 B | `9439a5fde12a3fb2` |
| `manuscript_figures/external/calibration_onestep.csv` | 1,295 B | `94a52725e255cd09` |
| `dripwater/HS4_dripwater_canonical.csv` | 61,286 B | `cc3de093b563a063` |
| `dr_app/HS4_example_inputs/` | (folder, 4 files) | `` |

Files to add that are NOT in the repository:

- `manuscript_figures/external/drip_rate_realisations_ap.csv.gz`, the age-propagated realisation ensemble (Supplementary Methods 7, Supplementary Fig. 12). Take it from the existing v1.0.0-beta deposit (Drip_Rate_Final.zip, 578 MB) or regenerate with Dr Paleo in age mode at the production settings (dr_app/PRODUCTION_SETTINGS.md rev 4).
- The DGT deployment table behind Supplementary Fig. 20 (Salmanzadeh et al., in preparation), released as Source Data with SF20.
- The Source Data workbook for the Supplementary Figures if it is kept separate from HS4_SourceData.xlsx.

Suggested metadata for the new version: title "PaleodripRates: data for the HS4 (Heshang Cave) drip-rate and precipitation reconstruction"; upload type Dataset; version 1.1.0; all ten authors of the article; licence CC BY 4.0 (as the existing version); related identifier "isSupplementedBy" the code release DOI and "isSupplementTo" the article DOI once known; description naming the resubmission and pointing to PATCH_NOTES.md in the code release.