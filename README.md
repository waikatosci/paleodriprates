# 💧 Dr Paleo — Paleodriprates

Kinetic proxy for stalagmite drip rate and precipitation reconstruction.

This repository contains code, data, and **Dr Paleo** — a browser-based application for reconstructing cave drip rates and Holocene precipitation from stalagmite trace metals. The kinetic proxy model exploits the dissociation kinetics of organic-metal complexes (OMCs) bound to transition metals (Co, Ni, Cu, Zn, Al, and others) in cave dripwater: because drip rate governs the thin-film residence time available for OMC dissociation, trace metal concentrations in stalagmite calcite directly encode past drip rates (drips min⁻¹). Through calibrated site-specific regressions and Monte Carlo propagation of uncertainties, these are translated into absolute precipitation estimates (mm yr⁻¹).

The methods are described in:

> Hartland, A., Goswami, B., Park, J., Höpker, S.N., Torres Rojas, D., Liao, J., Fox, B.R.S., Marwan, N., Breitenbach, S.F.M. & Hu, C. (submitted). Quantitative Holocene precipitation reconstruction from stalagmite trace metal kinetics reveals East Asian monsoon drivers. *Nature Geoscience*. DOI: [pending].

---

## Quick Start

See **[dr_app/QUICKSTART.md](dr_app/QUICKSTART.md)** for a step-by-step guide, including a worked example with the HS4 dataset.

```bash
git clone https://github.com/waikatosci/paleodriprates.git
cd paleodriprates/dr_app
pip install flask numpy pandas scipy
python app.py
# Open http://localhost:5000 → 💧 Dr Paleo
```

---

## Repository Structure

There are **two ways to run the model** — a Flask web app (recommended for most users) and a command-line driver (for batch processing). Both share a common core of forward-model code at the repo root.

```
paleodriprates/
│
├── README.md                          ← This file
├── LICENSE                            ← MIT License
├── requirements.txt
├── drip_rate.yml                      ← Conda environment file
│
├── # ---- Shared core (used by both the CLI and the web app) ----
├── model.py                           ← Forward model: h(V) trace element kinetics
├── params.py                          ← Physical constants (VMAX, VMIN, VRES, etc.)
├── utils.py                           ← Generic helpers (progress bars, IO)
├── drip_rate_util.py                  ← Drip rate–specific helpers (outlier
│                                          detection, residual optimisation)
│
├── # ---- CLI entry point ----
├── drip_rate_mc_realisations.py       ← Canonical CLI driver. Parallel Monte
│                                          Carlo over the kinetic inversion;
│                                          emits both percentile summaries and
│                                          the full MC realisation ensemble used
│                                          by the RQA stationarity tests.
│                                          (Reads Drip_rate.xlsx at the root.)
├── Drip_rate.xlsx                     ← Reference dataset: depth–age and trace
│                                          element profiles for CLI runs.
├── drip_rate_stationarity_tests.py    ← RQA stationarity tests on the MC
│                                          realisation ensemble emitted above.
│
├── dr_app/                            ← 💧 Dr Paleo Flask web application
│   ├── app.py                         ← Flask app entry point
│   ├── model_stochastic.py            ← Stochastic-prior wrapper around model.py
│   ├── driprates_stochastic.py        ← Drip rate PDF computation with priors
│   ├── concentration_prior.py         ← Log-normal priors for [TE]aq and [Ca]aq
│   ├── QUICKSTART.md                  ← Step-by-step guide for the web app
│   ├── launch_windows.bat             ← Windows one-click launcher
│   ├── launch_mac_linux.sh             ← macOS/Linux one-click launcher
│   ├── HS4_example_inputs/            ← Example data (Heshang Cave stalagmite HS4)
│   │   ├── HS4_age_depth.csv          ← U-Th dating table
│   │   ├── HS4_TE.csv                 ← Co and Ni LA-ICP-MS profiles (CSV)
│   │   └── HS4_TE.xlsx                ← Same data, Excel format
│   ├── uploads/                       ← Uploaded CSV data (created at runtime,
│   │                                      gitignored)
│   └── outputs/                       ← Model outputs and cached proxy records
│                                          (created at runtime, gitignored)
│
├── bayprox/                           ← BayProX: Bayesian proxy–age modelling
│   ├── data.py
│   ├── agedepth.py
│   ├── proxyrecord.py
│   ├── simulate.py
│   ├── visualize.py
│   └── other/
│       ├── motabar/                   ← Vendored MoTaBaR (Heitzig, PIK Potsdam):
│       │                                  Monotonic Tail-Adapted Bayesian
│       │                                  Regression. Live dependency of
│       │                                  bayprox/agedepth.py.
│       └── calibration/               ← Radiocarbon calibration data: IntCal09,
│                                          IntCal13, Marine13, Hua-Barbetti
│                                          postbomb tables.
│
├── companion_analysis/
│   └── RQA_HS4_ensemble.py            ← Recurrence quantification analysis on
│                                          the HS4 MC ensemble.
│
├── precip_recon/                      ← Precipitation reconstruction notebook
│   ├── P_quantification_Holocene.ipynb
│   ├── precip_recon_readme.txt
│   └── (regression input/output Excel files)
│
└── legacy/                            ← Read-only archive of superseded files
    ├── README.md                          (provenance + rationale)
    ├── scripts/                           (earlier CLI variants)
    ├── env/                               (Windows 7 conda env file)
    └── bayprox/                           (dated bayprox snapshots, scratch
                                            files)
```

### Two front-ends, one core

The **shared core** (`model.py`, `params.py`, `utils.py`, `drip_rate_util.py` at the repo root) contains the forward kinetic model and the helpers needed to run it. Both front-ends import from these directly:

- **Web app** (`dr_app/app.py`) — implements its own Monte Carlo loop on top of the shared core, with stochastic priors (via `model_stochastic.py`, `driprates_stochastic.py`, `concentration_prior.py`). Reads CSVs uploaded through the browser.
- **CLI driver** (`drip_rate_mc_realisations.py`) — parallel MC over the shared core, no stochastic priors. Reads `Drip_rate.xlsx` at the repo root.

The two front-ends emit comparable outputs (percentile summaries + full MC realisation ensemble), but they are **separate implementations** of the MC loop. If you change the kinetic model (`model.py`) you change both. If you change one front-end's MC handling, you do **not** automatically affect the other.

---

## Dependencies

**Python 3.9+** is required.

```
numpy
pandas
scipy
matplotlib
openpyxl
Pillow
progressbar2
flask                  # Dr Paleo web app
```

`bayprox` is a custom Bayesian proxy–age modelling library included in this repository. No external installation is required.

---

## Usage Pathways

### 1. 💧 Dr Paleo — Browser-based (recommended)

The web application provides a complete GUI for the reconstruction pipeline.

```bash
cd dr_app
python app.py
```

Open `http://localhost:5000` in any modern browser. Dr Paleo walks you through data upload, parameter configuration, model execution, and result visualisation. See [QUICKSTART.md](QUICKSTART.md) for a full walkthrough.

**Windows:** Double-click `launch_windows.bat`
**macOS/Linux:** Run `./launch_mac_linux.sh`

### 2. Command-line driver

For batch processing or integration into existing pipelines:

```bash
python drip_rate_mc_realisations.py
```

Reads input from `Drip_rate.xlsx` at the repo root. Emits drip-rate percentile summaries **and** the full Monte Carlo realisation ensemble (CSV) used by `drip_rate_stationarity_tests.py`. Parallel by default; uses Python's `concurrent.futures` thread pool.

> Earlier CLI variants (`Drip_rate.py`, `Drip_rate_serial.py`, `Drip_rate_parallel.py`) are preserved in [`legacy/scripts/`](legacy/) for provenance — they are no longer maintained. The current driver was previously named `Drip_rate_parallel_fr.py`; the `_fr` suffix stood for "full realisations".

### 3. Precipitation reconstruction

After obtaining drip rate percentiles (via Dr Paleo or command-line):

```bash
cd precip_recon
jupyter notebook P_quantification_Holocene.ipynb
```

Chains site-specific regressions (drip rate → discharge → precipitation) with Monte Carlo propagation.

### 4. RQA stationarity tests

For testing stationarity of the Monte Carlo realisation ensemble:

```bash
python drip_rate_stationarity_tests.py
```

Runs five tests (ADF, KPSS, Mann-Kendall, Ljung-Box, Runs) on the realisation CSV exported by Dr Paleo.

---

## 💧 Dr Paleo — Web Application

### Six-panel workflow

| Panel | Purpose |
|-------|---------|
| **Data Inputs** | Upload CSV files (depth/age, trace elements, isotopes, aqueous monitoring). Map columns and set units. |
| **Model Parameters** | Set cave conditions (temperature, Ca, drip rate), kinetic parameters (Kp, Kd, fractions), and aqueous chemistry per TE. |
| **Analysis Mode** | Choose full quantification (absolute drips min⁻¹) or semi-quantitative (% of reference). |
| **Output Options** | Configure V-grid, realisations, proxy record caching. |
| **Run** | Execute the model with progress bar, ETA, and live log. |
| **Results** | Interactive charts (time series, heatmap, Smart & Friedrich classification, age model) and downloadable CSVs. |

### Supported elements

Cu, Ni, Co (full theoretical Kp via Wang & Xu 2001 + Lindeman 2022 empirical), Zn, Cd, Pb, V, Mn, Fe, Al (literature Kp), and user-defined elements.

### Output files

| File | Contents |
|------|----------|
| `drip_rate_summary.csv` | Percentile summary (pc05–pc95) at each timestep |
| `drip_rate_realisations.csv` | Full MC ensemble for RQA analysis |
| `age_model.csv` | Depth–age mapping with errors |
| `input_summary.csv` | All input parameters for reproducibility |

---

## Statistical Tests

`drip_rate_stationarity_tests.py` performs five stationarity/trend tests on Monte Carlo realisations:

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

Raw proxy data and monitoring records are included in the Excel files.
Full archive (including all calibration datasets): [Zenodo DOI: 10.5281/zenodo.16392750](https://doi.org/10.5281/zenodo.16392750).

---

## Citation

**Manuscript:**
> Hartland, A., Goswami, B., Park, J., Höpker, S.N., Torres Rojas, D., Liao, J., Fox, B.R.S., Marwan, N., Breitenbach, S.F.M. & Hu, C. (2025). Quantitative Holocene precipitation reconstruction from stalagmite trace metal kinetics reveals East Asian monsoon drivers. *Nature Geoscience*. DOI: [pending].

**Repository:**
> Hartland, A. et al. (2025). PaleodripRates: Code for stalagmite drip rate and precipitation reconstruction. Zenodo. https://doi.org/10.5281/zenodo.16392750

**Partitioning model foundation:**
> Lindeman, I., Hansen, M., Scholz, D., Breitenbach, S.F.M. & Hartland, A. (2022). Effects of organic matter complexation on partitioning of transition metals into calcite: Cave-analogue crystal growth experiments. *Geochimica et Cosmochimica Acta* 317, 118–137.

**OMC conceptual basis:**
> Hartland, A. & Zitoun, R. (2018). Transition metal availability to speleothems controlled by organic binding ligands. *Geochemical Perspectives Letters* 8, 22–25.

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

## Acknowledgments

Funded by EU Horizon 2020 Marie Skłodowska-Curie Actions (no. 691037, QUEST — QUantitative paleoEnvironments from SpeleoThems), Te Apārangi Royal Society of New Zealand (RIS-UOW1501), Ministry for Business, Innovation and Employment (UOWX2102), and a Rutherford Discovery Fellowship (RDF-UOW1601).

For questions or contributions, open an issue on GitHub or contact the corresponding author: [waikatoscientific@gmail.com](mailto:waikatoscientific@gmail.com)
