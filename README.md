# PaleodripRates

Kinetic proxy for stalagmite drip rate and precipitation reconstruction.

This repository contains code, data, and a browser-based application for reconstructing cave drip rates and Holocene precipitation from stalagmite trace metals. The kinetic proxy model exploits the dissociation kinetics of organic-metal complexes (OMCs) bound to transition metals (Co, Ni) in cave dripwater: because drip rate governs the thin-film residence time available for OMC dissociation, trace metal concentrations in stalagmite calcite directly encode past drip rates (drips min⁻¹). Through calibrated site-specific regressions and Monte Carlo propagation of uncertainties, these are translated into absolute precipitation estimates (mm yr⁻¹).

The methods are described in: Hartland, A., Goswami, B., Park, J., Höpker, S.N., Torres Rojas, D., Liao, J., Fox, B.R.S., Marwan, N., Breitenbach, S.F.M. & Hu, C. (2025). Quantitative Holocene precipitation reconstruction from stalagmite trace metal kinetics reveals East Asian monsoon drivers. *Nature Geoscience*. DOI: [pending].

---

## Repository Structure

```
paleodriprates/
├── model.py                   # Core dissociation kinetics and expectation integrals
├── params.py                  # Global parameters (V grid, BayProX settings, outlier thresholds)
├── drip_rate_util.py          # Bayesian inversion for drip rates (scalar path)
├── Drip_rate_parallel.py      # Main driver: Excel → BayProX → drip rate PDFs (parallel)
├── Drip_rate_serial.py        # Serial alternative (slower, same outputs)
├── Drip_rate.py               # Legacy single-threaded driver
├── utils.py                   # Progress bars and helpers
├── Drip_rate.xlsx             # Input/output Excel workbook (HS4 data)
├── requirements.txt           # Python dependencies
├── LICENSE                    # MIT License
│
├── bayprox/                   # BayProX — Bayesian proxy record estimation library
│   ├── __init__.py
│   ├── agedepth.py            # Age-depth modelling
│   ├── data.py                # Data structures (SampleInfo)
│   ├── main.py                # Proxy record computation
│   ├── proxyrecord.py         # ProxyDistributions class
│   ├── simulate.py            # Simulation utilities
│   └── visualize.py           # Plotting routines
│
├── companion_analysis/
├── RQA_HS4_ensemble.py        # Windowed ensemble RQA
├── drip_rate_stationarity_tests.py    # Statistical tests for event significance and stationarity 
│
├── precip_recon/              # Precipitation reconstruction (chained regressions)
│   ├── P_quantification_Holocene.ipynb   # Jupyter notebook: Monte Carlo P reconstruction
│   ├── Heshang-Yichang_PT_regression_data.xlsx
│   ├── T_recon_Wang_et_al.xlsx
│   ├── drip_rate_percentiles.xlsx
│   ├── p_reconstruction.csv   # Output: annual P posteriors (median, pc25, pc75)
│   ├── p_plot.png             # Output: Holocene P figure
│   └── precip_recon_readme.txt
│
└── dr_app/                    # Browser-based Drip Rate Estimator (Flask)
    ├── app.py                 # Flask application (single-file, embedded HTML/JS)
    ├── concentration_prior.py # Lognormal priors for stochastic aqueous chemistry
    ├── driprates_stochastic.py# Drop-in driprates() with concentration uncertainty
    ├── model_stochastic.py    # Vectorised forward model (precomputed E₁/E₂)
    ├── uploads/               # User-uploaded CSV files (created at runtime)
    └── outputs/               # Model output files (created at runtime)
```

---

## Dependencies

Python 3.12+ with the following packages:

numpy, pandas, scipy, matplotlib, openpyxl, Pillow, progressbar

The `bayprox` library is included in the repository. The web application additionally requires Flask (`pip install flask`).

## Installation

```bash
git clone https://github.com/waikatosci/paleodriprates.git
cd paleodriprates
pip install -r requirements.txt
```

---

## Usage

### 1. Command-Line Workflow (Excel-Based)

The original workflow reads input data from `Drip_rate.xlsx` and writes results back to the same workbook.

**Input sheets:** `0.Settings` (station name, calibration range, plot options), `1.Depth_Age` (U-Th dated depth-age pairs), `2.Trace_Elems` (Co/Ni concentrations in ppm, kinetic parameters, aqueous chemistry in ppb), `3.Isotopes` (δ¹⁸O, δ¹³C).

**Run the model:**

```bash
python Drip_rate_parallel.py    # recommended (multi-threaded)
python Drip_rate_serial.py      # alternative (single-threaded)
```

**Output sheets:** `4.Output` (prior remainder variance, unified depth grid, proxy record statistics), `5.OutDripRate` (drip rate percentiles: pc05–pc95 and median at each year BP), `6.OutIsotope` (median δ¹⁸O vs age), `7.PlotDripRate` and `8.PlotDripIsotope` (embedded figures).

Typical runtime is 30–60 minutes depending on the number of trace elements and age range.

### 2. Web Application (dr_app)

The `dr_app/` folder provides a browser-based interface that replaces the Excel workflow. It accepts CSV uploads, exposes all model parameters through an interactive UI, and produces downloadable output files including full Monte Carlo realisation ensembles suitable for recurrence quantification analysis (RQA).

**Quick start:**

```bash
cd dr_app
pip install flask        # one-time
python app.py            # opens at http://localhost:5000
```

**Workflow:**

1. **Upload** CSV files for depth-age, trace element(s), and optionally isotopes.
2. **Configure** model parameters: dissociation rates (ln k_d, σ), partition coefficients (K_p), fast/inert fractions, aqueous concentrations, cave temperature.
3. **Run** the model with real-time progress tracking in the browser.
4. **Download** outputs: `drip_rate_summary.csv` (percentile time series), `drip_rate_realisations.csv` (full ensemble for RQA), `age_model.csv`, and `input_summary.csv` (parameter archive for reproducibility).

**Stochastic aqueous chemistry:** The web app supports lognormal concentration priors for both trace element and calcium aqueous concentrations, propagating solution chemistry uncertainty through the forward model. When enabled, the stochastic engine (`driprates_stochastic.py` + `model_stochastic.py`) precomputes the expensive E₁/E₂ integrals once over the V grid, then marginalises over N concentration samples by rescaling K₀ — reducing cost by approximately 100× compared with repeated full integration. This mode is configured via the "Concentration priors" panel in the UI.

**Cached proxy records:** If a previous run produced `ProxyRecord.pkl`, the app can reload it to skip the BayProX computation (~20 min), going directly to drip rate estimation.

### 3. Precipitation Reconstruction

After obtaining drip rate percentiles (from either workflow above), run the Jupyter notebook to compute absolute annual precipitation:

```bash
cd precip_recon
jupyter notebook P_quantification_Holocene.ipynb
```

This chains site-specific regressions (Heshang drip rate → P/T → Yichang P/T → P) with Monte Carlo sampling (N = 10,000 per age horizon) of drip rate posteriors, regression residuals, and temperature uncertainty (Wang et al. 2018; RMSE = 2.6°C). Outputs are `p_reconstruction.csv` (median, pc25, pc75) and `p_plot.png`.

### 4. Recurrence Quantification Analysis

The 1,000-realisation ensemble CSVs exported by the web application (or generated by a custom extraction script) serve as input for windowed ensemble RQA. The RQA analysis code is maintained separately (see `RQA_HS4_ensemble.py` in the companion analysis repository). Parameters used in the published analysis: τ = 1, m = 5 (estimated globally via AMI and Cao's method), RR = 0.05, Theiler window = 4, sliding window = 100 points (step 5). DET and TRANS are reported with 5th–95th percentile uncertainty envelopes from the full ensemble.

---

## Key Modules

**model.py** — Core kinetic model. Computes the expectation integrals E₁(V) and E₂(V) over a log-normal distribution of OMC dissociation rate constants, the forward function h(V) mapping drip rate to calcite trace metal concentration, and the Jacobian h′(V) for PDF transformation. Also provides `kp_theory()` for theoretical partition coefficients (Wang & Xu 2001) and `dr_pdfseries()` for the full drip rate PDF time series.

**drip_rate_util.py** — Bayesian inversion wrapper. Calls `model.te_pdfseries()` to obtain the BayProX proxy PDFs, converts aqueous concentrations from ppb to mol/L, and passes everything to `model.dr_pdfseries()`. Returns V_pdf (drip rate PDF matrix), age axis, and V_span.

**concentration_prior.py** — Lognormal priors for aqueous concentrations. Can be constructed from a monitoring time series (`from_series`), manual mean/SD (`from_mean_sd`), or direct lognormal parameters. Sampling always returns ppb. Used by the stochastic engine to propagate solution chemistry uncertainty.

**driprates_stochastic.py** — Drop-in replacement for `drip_rate_util.driprates()`. When the TE dictionary contains `ca_prior` and/or `aq_prior` keys (ConcentrationPrior instances), it samples N concentration pairs and marginalises via geometric-mean log-PDF averaging. Otherwise delegates to the original scalar path.

**model_stochastic.py** — Vectorised forward model for stochastic mode. Precomputes E₁(V) and E₂(V) once on the full V grid, then for each concentration sample rescales h(V) = K₀ × (1 − n_S × E₁) without re-integrating — yielding ~100× speedup over naive repeated calls.

---

## Data

Input data for the HS4 stalagmite (Heshang Cave, China) are provided in `Drip_rate.xlsx`. Trace element concentrations in the carbonate are in ppm (solid phase); aqueous monitoring concentrations are in ppb. The full dataset — including raw ICP-MS data, monitoring records, and all model outputs — is archived on Zenodo: [DOI: 10.5281/zenodo.16392750](https://doi.org/10.5281/zenodo.16392750).

## Citation

If using this code or data, please cite:

Hartland, A., Goswami, B., Park, J., Höpker, S.N., Torres Rojas, D., Liao, J., Fox, B.R.S., Marwan, N., Breitenbach, S.F.M. & Hu, C. (2025). Quantitative Holocene precipitation reconstruction from stalagmite trace metal kinetics reveals East Asian monsoon drivers. *Nature Geoscience*. DOI: [pending].

For the repository: Hartland, A. et al. (2025). PaleodripRates: Code for stalagmite drip rate and precipitation reconstruction. Zenodo. https://doi.org/10.5281/zenodo.16392750

## License

MIT License — see [LICENSE](LICENSE) for details.

## Acknowledgements

Funded by EU Horizon 2020 Marie Skłodowska-Curie Actions (grant no. 691037; QUEST); Te Apārangi Royal Society of New Zealand (RIS-UOW1501); Ministry for Business, Innovation and Employment (UOWX2102); and a Rutherford Discovery Fellowship (RDF-UOW1601).

For questions or contributions, open an issue or contact the corresponding author (adam.hartland@waikato.ac.nz).
