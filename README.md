# Dr Paleo: Paleodriprates

Kinetic proxy for stalagmite drip rate and precipitation reconstruction.

This repository contains the code, data and Dr Paleo, a browser-based application for reconstructing cave drip rates and Holocene precipitation from trace metals in stalagmites. The kinetic proxy rests on the dissociation of organic-metal complexes (OMCs) that carry transition metals (Co, Ni, Cu, Zn and others) in cave dripwater. Drip rate sets how long the water film on the stalagmite surface has for these complexes to dissociate, so the metal content of the calcite records past drip rate (drips min⁻¹). Calibrated site-specific regressions, with uncertainties propagated by Monte Carlo, convert drip rate to absolute precipitation (mm yr⁻¹).

The methods and the Heshang Cave (HS4) record are described in:

> Hartland, A., Goswami, B., Park, J., Höpker, S.N., Torres Rojas, D., Liao, J., Fox, B.R.S., Marwan, N., Breitenbach, S.F.M. & Hu, C. (submitted 2026, Nature Communications; NCOMMS-26-041445-T). Decoupled infiltration and isotope signals reveal a hidden East Asian monsoon megadrought. Preprint and DOI to follow.

---

## Quick Start

See [dr_app/QUICKSTART.md](dr_app/QUICKSTART.md) for a step-by-step guide with a worked HS4 example, and [dr_app/PRODUCTION_SETTINGS.md](dr_app/PRODUCTION_SETTINGS.md) for the settings that reproduce the published record.

```bash
git clone https://github.com/waikatosci/paleodriprates.git
cd paleodriprates/dr_app
pip install flask numpy pandas scipy
python app.py
# Open http://localhost:5000
```

---

## Repository Structure

There are two ways to run the model: a Flask web app (recommended for most users) and a command-line driver (for batch processing). Both use the same forward-model code at the repository root. The remaining folders hold the calibration, the analyses behind the paper, and the scripts that draw every figure.

```
paleodriprates/
│
├── README.md                          This file
├── PATCH_NOTES.md                     Change log for the revision (what changed and why)
├── LICENSE                            MIT License
├── requirements.txt                   pip install path
├── drip_rate.yml                      Conda environment file (recommended)
│
│   ---- Shared core (used by both the CLI and the web app) ----
├── model.py                           Forward model: h(V) trace-element kinetics
├── params.py                          Physical constants (VMAX, VMIN, VRES, etc.)
├── utils.py                           Generic helpers (progress bars, IO)
├── drip_rate_util.py                  Drip-rate helpers (outlier detection,
│                                          residual optimisation)
│
│   ---- CLI entry point ----
├── drip_rate_mc_realisations.py       Canonical CLI driver: parallel Monte Carlo
│                                          over the kinetic inversion; writes
│                                          percentile summaries and the full
│                                          realisation ensemble. Reads Drip_rate.xlsx.
├── Drip_rate.xlsx                     Reference dataset: age-depth, trace elements
│                                          and isotopes for CLI runs
├── drip_rate_stationarity_tests.py    Stationarity tests on the realisation ensemble
│
├── dr_app/                            Dr Paleo Flask web application
│   ├── app.py                         Flask app entry point
│   ├── model_stochastic.py            Stochastic-prior wrapper around model.py
│   ├── driprates_stochastic.py        Drip-rate PDF computation with priors
│   ├── concentration_prior.py         Log-normal priors for [TE]aq and [Ca]aq
│   ├── QUICKSTART.md                  Step-by-step guide for the web app
│   ├── PRODUCTION_SETTINGS.md         Settings that reproduce the HS4 record
│   ├── launch_windows.bat             Windows one-click launcher
│   ├── launch_mac_linux.sh            macOS/Linux one-click launcher
│   ├── HS4_example_inputs/            Heshang Cave stalagmite HS4
│   │   ├── HS4_age_depth.csv          ²³⁰Th ages (30) and the growth surface
│   │   ├── HS4_TE.csv / .xlsx         Full Ni and Co table as measured
│   │   │                                  (589 rows, both laboratories)
│   │   └── HS4_TE_canonical.csv       Reconstruction input: the 568 samples of
│   │                                      the primary (2009) run
│   ├── uploads/, outputs/             Created at runtime (gitignored)
│
├── calibration/                       Kd calibration against the monitoring record,
│                                          canonical-input builder, exclusion list,
│                                          below-resolution censoring, Fig. 3 curves
├── dripwater/                         HS4 dripwater monitoring record (2007 onward):
│                                          drip rate and full solution chemistry
│
├── companion_analysis/                Analyses cited in the Supplementary Information
│   ├── RQA_HS4_ensemble.py            Recurrence quantification of the MC ensemble
│   ├── run_single_metal.py            Ni-only and Co-only inversions
│   ├── ni_co_agreement.py             Ni-Co agreement and the 5.2 ka interval
│   ├── crosslab_*.py, run_crosslab_585.py
│   │                                  Second-laboratory (2019) run: matrix check,
│   │                                      sensitivity, Ca-corrected replication
│   ├── detrital_ternary_screen.py     Detrital and source-covariance screens
│   ├── events_differ_lithogenic.py,
│   │   make_fig_events_differ.py      8.2 ka and 5.2 ka event comparison
│   ├── dripwater_source_stability.py,
│   │   source_variation_propagation.py
│   │                                  Stability of the dripwater metal source
│   ├── element_covariance.py          Which elements covary with Ni and Co
│   └── dripwater_integration_timescale.py
│                                      Over what timescale the dripwater metals track drip
│                                          rate (Supplementary Methods 16)
│
├── extended_data/                     Scripts for Supplementary Figures 11-16
│                                          (see extended_data/README.md)
├── supplementary_figures/             Scripts for Supplementary Figures 1-5, and a
│                                          table mapping every Supplementary Figure
│                                          to the script that draws it
│
├── manuscript_figures/                Main-text figures and Source Data
│   ├── HS4_SourceData.xlsx            Single source for every display item
│   ├── generate_figures.py            Renders all main-text figures
│   ├── build_precip_onestep.py        One-step precipitation transfer
│   ├── ngeo_style.py                  Shared figure style
│   ├── figures/                       One script per figure (Figs 2-7)
│   ├── external/                      Bulk arrays and cited inputs (see its README)
│   └── output/                        Rendered figures and tables
│
├── precip_recon/                      Precipitation reconstruction notebook
│   ├── P_quantification_Holocene.ipynb
│   ├── precip_recon_readme.txt
│   └── (regression inputs and outputs)
│
├── bayprox/                           BayProX: Bayesian proxy-age modelling
│   ├── data.py, agedepth.py, proxyrecord.py, simulate.py, visualize.py
│   └── other/
│       ├── motabar/                   Vendored MoTaBaR (Heitzig, PIK Potsdam):
│       │                                  Monotonic Tail-Adapted Bayesian
│       │                                  Regression, used by bayprox/agedepth.py
│       └── calibration/               Radiocarbon calibration data (IntCal09,
│                                          IntCal13, Marine13, Hua-Barbetti)
│
└── legacy/                            Read-only archive of superseded files
    ├── README.md                      Provenance and rationale
    ├── scripts/                       Earlier CLI variants
    ├── env/                           Windows 7 conda environment file
    └── bayprox/                       Dated bayprox snapshots
```

### Two front-ends, one core

The shared core (`model.py`, `params.py`, `utils.py`, `drip_rate_util.py`) contains the forward kinetic model and the helpers needed to run it. Both front-ends import from it directly:

- Web app (`dr_app/app.py`): its own Monte Carlo loop on the shared core, with stochastic priors (`model_stochastic.py`, `driprates_stochastic.py`, `concentration_prior.py`). Reads CSV files uploaded through the browser.
- CLI driver (`drip_rate_mc_realisations.py`): parallel Monte Carlo on the shared core, without stochastic priors. Reads `Drip_rate.xlsx`.

The two front-ends write comparable outputs (percentile summaries and the full realisation ensemble) but are separate implementations of the Monte Carlo loop. A change to `model.py` affects both; a change to one front-end's Monte Carlo handling does not affect the other.

### HS4 inputs

The published record is the inversion of `HS4_TE_canonical.csv`: the 568 samples of the primary ICP-MS run (2009). The full table, `HS4_TE.csv`, also holds the 0.06 cm surface-cap point and the 20 samples of a second-laboratory run (2019), which are withdrawn from the reconstruction and examined separately (Supplementary Methods 14.3–14.4; `companion_analysis/crosslab_*.py`). `calibration/make_canonical_te_input.py` builds the canonical file from the full table and `calibration/excluded_points.csv`.

---

## Dependencies

Python 3.9 or later; the conda environment (`drip_rate.yml`) pins Python 3.12.

```
numpy
pandas
scipy
matplotlib
openpyxl
Pillow
progressbar
flask                  # Dr Paleo web app
```

Some companion analyses and figure scripts need further packages that the app does not use: `astropy` (spectral analysis, Supplementary Figure 15) and `cartopy` and `xarray` (site map, Fig. 2).

`bayprox` is a Bayesian proxy-age modelling library included in this repository; no separate installation is needed.

---

## Usage Pathways

### 1. Dr Paleo, browser-based (recommended)

The web application provides a complete interface for the reconstruction.

```bash
cd dr_app
python app.py
```

Open `http://localhost:5000` in any modern browser. Dr Paleo takes you through data upload, parameter settings, model runs and results. See [QUICKSTART.md](dr_app/QUICKSTART.md) for a full walkthrough.

- Windows: double-click `launch_windows.bat`
- macOS/Linux: run `./launch_mac_linux.sh`

### 2. Command-line driver

For batch processing or use in existing pipelines:

```bash
python drip_rate_mc_realisations.py
```

Reads `Drip_rate.xlsx` at the repository root and writes drip-rate percentile summaries and the full Monte Carlo realisation ensemble (CSV) used by `drip_rate_stationarity_tests.py`. Runs in parallel by default (Python's `concurrent.futures` thread pool).

> Earlier CLI variants (`Drip_rate.py`, `Drip_rate_serial.py`, `Drip_rate_parallel.py`) are kept in [`legacy/scripts/`](legacy/) for provenance and are no longer maintained. The current driver was previously named `Drip_rate_parallel_fr.py` (`_fr` for "full realisations").

### 3. Precipitation reconstruction

After obtaining drip-rate percentiles (from Dr Paleo or the command line):

```bash
cd precip_recon
jupyter notebook P_quantification_Holocene.ipynb
```

Chains the site-specific regressions (drip rate, discharge, precipitation) with Monte Carlo propagation.

### 4. Stationarity tests

```bash
python drip_rate_stationarity_tests.py
```

Runs five tests (ADF, KPSS, Mann-Kendall, Ljung-Box, runs) on a realisation ensemble written by Dr Paleo or the CLI driver.

### 5. Figures and companion analyses

```bash
cd manuscript_figures
python generate_figures.py          # main-text figures from HS4_SourceData.xlsx
```

Each Supplementary Figure is drawn by one script; `supplementary_figures/README.md` lists which. The companion analyses in `companion_analysis/` run from the repository root and write to `manuscript_figures/output/` and `manuscript_figures/external/`.

---

## Dr Paleo web application

### Six-panel workflow

| Panel | Purpose |
|-------|---------|
| Data Inputs | Upload CSV files (age-depth, trace elements, isotopes, dripwater monitoring). Map columns and set units. |
| Model Parameters | Set cave conditions (temperature, Ca, drip rate), kinetic parameters (Kp, Kd, fractions) and dripwater chemistry for each element. |
| Analysis Mode | Full quantification (absolute drips min⁻¹) or semi-quantitative (% of a reference). |
| Output Options | Configure the V-grid, number of realisations and proxy-record caching. |
| Run | Run the model with progress bar, time remaining and live log. |
| Results | Interactive charts (time series, heatmap, Smart & Friedrich classification, age model) and downloadable CSV files. |

### Supported elements

Cu, Ni, Co (theoretical Kp from Wang & Xu 2001 and empirical Kp from Lindeman et al. 2022), Zn, Cd, Pb, V, Mn, Fe, Al (literature Kp) and user-defined elements.

### Output files

| File | Contents |
|------|----------|
| `drip_rate_summary.csv` | Percentile summary (pc05-pc95) at each step |
| `drip_rate_realisations.csv` | Full Monte Carlo ensemble for RQA and stationarity tests |
| `age_model.csv` | Depth-age mapping with errors |
| `input_summary.csv` | All input parameters, for reproducibility |

---

## Statistical Tests

`drip_rate_stationarity_tests.py` performs five stationarity and trend tests on the Monte Carlo realisations:

| Test | Null hypothesis | Measures |
|------|----------------|----------|
| ADF | Unit root (non-stationary) | Mean stationarity |
| KPSS | Trend stationary | Trend stationarity |
| Mann-Kendall | No monotonic trend | Trend direction |
| Ljung-Box | No autocorrelation | Serial dependence |
| Runs | Random sequence | Non-randomness |

All tests report effect sizes and bootstrap 95% confidence intervals.

---

## Data Availability

The stalagmite proxy data, the U-Th chronology and the dripwater monitoring record are in this repository (`dr_app/HS4_example_inputs/`, `Drip_rate.xlsx`, `dripwater/`). Every main-text display item is drawn from the single workbook `manuscript_figures/HS4_SourceData.xlsx` (drip-rate reconstruction at σ = π/√6). The full trace-element measurement report of the second-laboratory run is in `manuscript_figures/external/HS4_TE_full_suite_both_labs.xlsx`. Bulk arrays too large for the repository (for example the age-propagated realisation ensemble) are in the Zenodo archive: [DOI 10.5281/zenodo.16392750](https://doi.org/10.5281/zenodo.16392750).

---

## Citation

Manuscript:
> Hartland, A., Goswami, B., Park, J., Höpker, S.N., Torres Rojas, D., Liao, J., Fox, B.R.S., Marwan, N., Breitenbach, S.F.M. & Hu, C. (submitted 2026, Nature Communications; NCOMMS-26-041445-T). Decoupled infiltration and isotope signals reveal a hidden East Asian monsoon megadrought. Preprint and DOI to follow.

Repository:
> Hartland, A. et al. (2025). PaleodripRates: Code for stalagmite drip rate and precipitation reconstruction. Zenodo. https://doi.org/10.5281/zenodo.16392750

Partitioning model foundation:
> Lindeman, I., Hansen, M., Scholz, D., Breitenbach, S.F.M. & Hartland, A. (2022). Effects of organic matter complexation on partitioning of transition metals into calcite: cave-analogue crystal growth experiments. Geochimica et Cosmochimica Acta 317, 118-137.

OMC conceptual basis:
> Hartland, A. & Zitoun, R. (2018). Transition metal availability to speleothems controlled by organic binding ligands. Geochemical Perspectives Letters 8, 22-25.

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---

## Acknowledgments

Funded by EU Horizon 2020 Marie Skłodowska-Curie Actions (no. 691037, QUEST, QUantitative paleoEnvironments from SpeleoThems), Te Apārangi Royal Society of New Zealand (RIS-UOW1501), the Ministry of Business, Innovation and Employment (UOWX2102) and a Rutherford Discovery Fellowship (RDF-UOW1601). The ICP-MS and ICP-AES analyses were funded by the National Natural Science Foundation of China (41731177) to C. Hu.

For questions or contributions, open an issue on GitHub or contact the corresponding author: [adam.hartland@lincolnagritech.co.nz](mailto:adam.hartland@lincolnagritech.co.nz)

---

## AI Assistance Statement

The scientific method, kinetic model and interpretations in this repository are the work of the authors. Generative AI (Anthropic's Claude) was used under the authors' direction for software engineering and documentation: refactoring and organising the codebase, building the Dr Paleo web interface, writing analysis and figure scripts for the revision, and checking that values, settings and figure provenance agree across the code, data and manuscript. All AI-assisted output was reviewed, tested and verified by the authors, who take full responsibility for the content and correctness of the code, data and documentation.
