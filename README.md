# PaleodripRates

**Kinetic proxy reconstruction of cave drip rates and Holocene precipitation from stalagmite trace metals**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16392750.svg)](https://doi.org/10.5281/zenodo.16392750)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## Overview

This repository contains the full Python codebase, data, and a browser-based web application for reconstructing past cave drip rates — and by extension Holocene precipitation — from stalagmite trace metal concentrations (Co, Ni).

The approach exploits the dissociation kinetics of transition metal–organic complexes (OMC) in cave dripwater. Because the rate at which metals dissociate from dissolved organic ligands is controlled by the residence time of water on the stalagmite surface (i.e., 1/drip rate), the concentration of Co and Ni incorporated into speleothem calcite carries a quantitative signal of past drip rate. Drip rates are inverted probabilistically using Bayesian age-depth modelling (BayProX), then chained to precipitation/temperature regressions and independent temperature proxies for fully uncertainty-bounded precipitation estimates (mm yr⁻¹).

The method is validated and applied to stalagmite HS4 from Heshang Cave, Yangtze region, China, reconstructing Holocene rainfall over the past ~9,500 years and revealing new insights into East Asian Summer Monsoon dynamics.

---

## Citation

If you use this code or data, please cite:

> Hartland, A., Goswami, B., Park, J., Höpker, S., Torres Rojas, D., Liao, J., Fox, B.R.S., Marwan, N., Breitenbach, S.F.M., & Hu, C. (submitted). **Quantitative Holocene precipitation reconstruction from stalagmite trace metal kinetics reveals East Asian monsoon drivers.** *Nature Geoscience.*

For the repository itself:

> Hartland, A. et al. (2025). *PaleodripRates: Code for stalagmite drip rate and precipitation reconstruction.* Zenodo. https://doi.org/10.5281/zenodo.16392750

The experimental foundation for the partitioning model (inorganic K_d values for Co, Ni, Cu) is:

> Lindeman, I., Hansen, M., Scholz, D., Breitenbach, S.F.M., & Hartland, A. (2022). Effects of organic matter complexation on partitioning of transition metals into calcite: Cave-analogue crystal growth experiments. *Geochimica et Cosmochimica Acta* 317, 118–137. https://doi.org/10.1016/j.gca.2021.10.032

The experimental foundation for the dissociation model (inert fractions for Ni and Co) is: 

> Sebastian N. Hopker, Sebastian F.M. Breitenbach, Megan Grainger, Claudine Stirling, Adam Hartland. (2024). Characterising the decay of organic metal complexes in speleothem-forming cave waters. *Geochimica et Cosmochimica Acta* 373, 98-108. https://doi.org/10.1016/j.gca.2024.03.024

The conceptual basis for OMC-controlled transition metal availability to speleothems is established in:

> Hartland, A. & Zitoun, R. (2018). Transition metal availability to speleothems controlled by organic binding ligands. *Geochemical Perspectives Letters* 8, 22–25. https://doi.org/10.7185/geochemlet.1826

---

## Repository Structure

```
paleodriprates/
│
├── model.py                        # Core dissociation kinetics and dr_pdfseries
├── drip_rate_util.py               # Bayesian inversion utilities, outlier detection
├── drip_rate_parallel.py           # Parallelised drip rate pipeline (recommended)
├── drip_rate_serial.py             # Serial version (for debugging / single-core)
├── params.py                       # Shared model parameters and constants
├── utils.py                        # Progress bar helpers
│
├── dr_app/
│   ├── app.py                      # Flask web application (self-contained)
│   ├── launch_windows.bat          # Windows one-click launcher
│   ├── launch_mac_linux.sh         # macOS/Linux one-click launcher
│   ├── uploads/                    # Uploaded CSV data (created at runtime)
│   └── outputs/                    # Model outputs and cached proxy records
│
├── bayprox/                        # BayProX: Bayesian proxy–age modelling library
│   ├── data.py
│   ├── agedepth.py
│   └── proxyrecord.py
│
├── data/
│   ├── Drip_rate.xlsx              # Reference dataset: depth–age, Co/Ni proxy records
│   └── Precip_from_drip_rates.xlsx # Precipitation regression inputs
│
├── precip_recon/
│   ├── P_quantification_Holocene.ipynb   # Precipitation reconstruction notebook
│   └── readme.txt                        # Supplementary notes
│
├── drip_rate_stationarity_tests.py       # Statistical tests for event significance and stationarity
├── RQA_HS4_ensemble.py                   # Ensemble RQA for drip rate and δ¹⁸O (Fig. 5)
├── figures/                        # Generated plots (sensitivity, reconstructions)
├── requirements.txt
└── LICENSE
```

---

## Dependencies

**Python 3.12+** is required.

```
numpy
pandas
scipy
matplotlib
openpyxl
Pillow
progressbar2
flask                  # web app only
```

`bayprox` is a custom Bayesian proxy–age modelling library included in this repository. No external installation is required.

Install all dependencies:

```bash
pip install -r requirements.txt
```

Or individually:

```bash
pip install numpy pandas matplotlib scipy openpyxl pillow progressbar2 flask
```

---

## Usage

There are two ways to run the model: via the **web application** (recommended for new users) or directly via the **Python scripts** (recommended for batch processing and customisation).

---

### Option 1 — Web Application

The web app provides a complete browser-based interface to the drip rate reconstruction pipeline. No terminal interaction is required after launch.

#### Launch

**Windows:** Double-click `dr_app/launch_windows.bat`

**macOS / Linux:**
```bash
cd dr_app
chmod +x launch_mac_linux.sh
./launch_mac_linux.sh
```

Then open your browser at **http://localhost:5000**

#### Workflow

The app guides you through six panels:

| Panel | Description |
|-------|-------------|
| **Data Input** | Upload CSV files for depth–age model, trace element 1 (e.g. Ni), trace element 2 (e.g. Co), and isotope record. Skip if using a cached proxy record. |
| **Model Parameters** | Configure geochemical parameters: element selection (Co, Ni, Cu, and others), molecular weight, Kp (defaults to speleothem-specific values from Lindeman et al. 2022, or −1 for theoretical calculation via Wang & Xu 2001), mean and std dev of ln(K_d), aqueous concentration, and solution fractions X_F (fast), X_I (inert), X_L (labile, auto-calculated). Cave temperature and Ca concentration are also set here. |
| **Output Options** | Set number of Monte Carlo realisations (default 1000) and enable/disable the cached proxy record. |
| **Run Model** | Execute the pipeline with real-time progress, ETA, and scrolling log output. |
| **Results** | Interactive chart showing the drip rate reconstruction with 5th–95th percentile uncertainty envelope. |
| **Downloads** | Download `drip_rate_summary.csv` (percentiles) and `drip_rate_realisations.csv` (full ensemble for RQA and other nonlinear analyses). |

#### Cached proxy record

The most computationally expensive step is the BayProX proxy–age modelling (~20 min). After the first run, `ProxyRecord.pkl` is saved to `dr_app/outputs/`. On subsequent runs, tick **"Use cached proxy record"** to skip directly to the drip rate calculation (completes in minutes).

To use a proxy record from a prior script-based run, copy `ProxyRecord.pkl` from your working directory into `dr_app/outputs/` before launching the app.

---

### Option 2 — Python Scripts

#### Drip rate reconstruction

Run the parallelised version for best performance:

```bash
python drip_rate_parallel.py
```

Or the serial version:

```bash
python drip_rate_serial.py
```

**Inputs** (configured in `params.py` and read from `Drip_rate.xlsx`):
- Depth–age control points (U/Th dates) for BayProX age modelling
- Trace element proxy records (Co/Ca, Ni/Ca) as depth series
- Stable isotope record (δ¹⁸O) as depth series

**Outputs** (written to `Drip_rate.xlsx`):
- Drip rate PDFs (median + percentile envelopes) — `PlotDripRate` sheet
- Isotope comparison — `PlotDripIsotope` sheet
- `ProxyRecord.pkl` — cached BayProX output for reuse
- `drip_rate_realisations.csv` — full Monte Carlo ensemble (N_timesteps × N_realisations)

#### Precipitation reconstruction

After obtaining drip rate outputs, open the Jupyter notebook:

```bash
jupyter notebook precip_recon/P_quantification_Holocene.ipynb
```

This performs chained Monte Carlo regression linking drip rate → precipitation, propagating uncertainty from drip rate posteriors, monitoring regressions, and independent temperature estimates. Outputs include `p_reconstruction.csv` (annual precipitation posteriors) and `p_plot.png`.

---

## Model Description

### Kinetic proxy — organic metal complex dissociation

Transition metals (Co²⁺, Ni²⁺) in cave dripwater are predominantly bound to dissolved natural organic matter (NOM) as metal–ligand complexes (ML). The rate at which free M²⁺ becomes available for incorporation into calcite is governed by first-order dissociation kinetics:

```
[ML](t) = [ML]₀ · exp(−k · t)
```

where `k` is the dissociation rate constant and `t` is the thin-film residence time on the stalagmite surface, which is approximately 1/drip rate. Slower drip rates → longer residence → more dissociation → higher M/Ca in calcite. This provides the quantitative link between stalagmite trace metal concentrations and past drip rate.

### Solution fractions

Total aqueous metal is partitioned into:
- **X_I** — inert fraction: not bioavailable, does not participate in OMC dissociation
- **X_F** — fast fraction: labile but dissociates too rapidly to fractionate with drip rate
- **X_L** — labile fraction: undergoes drip-rate-sensitive OMC dissociation (X_L = 1 − X_I − X_F)

### Partition coefficients

Speleothem-specific inorganic K_d values (Lindeman et al. *GCA* 2022, cave-analogue GeoMic experiments):

| Element | K_d (inorganic) | K_d (+ NOM) | Notes |
|---------|----------------|-------------|-------|
| Co | 4.4 | 0.41 | PCP-sensitive; NOM strongly suppresses incorporation |
| Ni | 1.1 | 0.029 | PCP-insensitive (K_d ≈ 1); **recommended primary proxy** |
| Cu | 44 | 0.92 | Very high inorganic K_d; use with caution |

Ni is the preferred proxy element because its K_d ≈ 1 means prior calcite precipitation (PCP) has negligible effect on Ni/Ca, making it insensitive to aquifer drying artefacts that affect alkaline earth metal proxies (Sr, Mg).

### Theoretical K_p (Wang & Xu 2001)

When Kp = −1, the model calculates the partition coefficient from cave temperature using the lattice strain model:

```
log Kp = [ a(ΔGn_M − ΔGn_Ca) + b(r_M − r_Ca) − (ΔGf_M − ΔGf_Ca) ] / (−2.303·R·T)
```

For Co and Ni, the empirical values from Lindeman et al. (2022) are used in preference to the theoretical calculation, as the lattice strain model overestimates K_d for these elements under speleothem conditions.

---

## Statistical Tests

`drip_rate_stationarity_tests.py` runs five statistical tests against the Monte Carlo realisation ensemble and the full-posterior summary. Set `CSV_PATH` and `XLSX_PATH` at the top of the script to point to your drip rate outputs, then:

```bash
python drip_rate_stationarity_tests.py
```

**Inputs:**
- `drip_rate_realisations.csv` — full Monte Carlo ensemble (age column + one column per realisation)
- `Drip_rate_data_HS4.xlsx` — full-posterior summary with age, median, p25, p75 columns

**Tests performed:**

| # | Test | Windows | Method |
|---|------|---------|--------|
| 1 | **5.2 ka aridity pulse** | Event 4,500–5,200 BP vs. background 5,500–6,500 BP | IQR non-overlap; pooled MWU + KS; per-realisation MWU |
| 2 | **8.2 ka anomaly** | Event 8,000–8,500 BP vs. pre-event 9,000–9,500 BP and post-event 7,500–8,000 BP | IQR separation; pooled MWU + KS; per-realisation MWU |
| 3 | **Post-1820 CE surge** | Modern −400 to 50 BP vs. late-Holocene baseline 500–2,000 BP | Pooled MWU + KS; IQR overlap |
| 4 | **Holocene trend** | Full series | Mann–Kendall τ + Theil–Sen slope; early vs. late MWU |
| 5 | **Stationarity** | Full series | Augmented Dickey–Fuller (ADF) unit-root test |

All two-sample tests report rank-biserial effect size and bootstrap 95% CIs on the median difference (10,000 iterations; subsampled to n = 50,000 per group for large pools). Per-realisation MWU loops print progress every 100 realisations. ADF requires `statsmodels` (`pip install statsmodels`).

---

## Recurrence Quantification Analysis (RQA)

`RQA_HS4_ensemble.py` computes windowed ensemble RQA for both the HS4 drip rate and δ¹⁸O Monte Carlo ensembles, producing Fig. 5 and the supplementary recurrence plot figure (FigS1).

```bash
python RQA_HS4_ensemble.py
```

**Inputs:**

| File | Description |
|------|-------------|
| `drip_rate_realisations.csv` | Age column + r0…r999 drip rate realisations (drips min⁻¹) |
| `d18O_realisations.csv` | Age column + r0…r999 δ¹⁸O realisations (‰ VPDB) |
| `Drip_rate.xlsx` sheet `5.OutDripRate` | Drip rate summary: age, pc05–pc95, median |
| `Drip_rate.xlsx` sheet `6.OutIsotope` | δ¹⁸O summary: age, median |

**Outputs:**

| File | Description |
|------|-------------|
| `rqa_ensemble_results.csv` | Per-window DET and TRANS: median, 5th and 95th percentile for both proxies |
| `rqa_parameters.csv` | Global embedding parameters used in the run (τ, m, Theiler window, RR, l_min, window size/step) |
| `Fig5_RQA_HS4.pdf/.png/.eps` | Main figure — four stacked panels (Nature Geoscience style, 8.5 × 7 in) |
| `FigS1_RP_HS4.pdf/.png` | Supplementary full-record recurrence plots for δ¹⁸O and drip rate |

**Figure layout (Fig. 5):**
- **Panel A** — RQA metrics (DET, TRANS with 5–95th percentile envelopes) for δ¹⁸O
- **Panel B** — δ¹⁸O ensemble PDF heatmap with median overlay (cividis, inverted y-axis)
- **Panel C** — RQA metrics for drip rate
- **Panel D** — Drip rate ensemble PDF heatmap with median and IQR overlays

Dashed vertical markers are drawn automatically on panels A/C and B/D wherever DET falls and stays below `DET_THRESHOLD` (default 0.4) for at least `MIN_SUSTAIN` (default 80) consecutive windows, flagging sustained dynamical regime transitions.

**Method:**

Two RQA measures are computed in a sliding window (100 points, step 5) on each of the 1,000 realisations independently:

- **DET (Determinism)** — fraction of recurrent points forming diagonal lines of length ≥ 2 (Zbilut & Webber 1992). Higher DET indicates more periodic, predictable dynamics.
- **TRANS (Transitivity)** — recurrence network clustering coefficient (Donner et al. 2010). Reflects the geometric density of the reconstructed attractor; high values indicate constrained, low-dimensional dynamics.

**Embedding parameters** (τ and m) are estimated *once* from the full median series of each proxy, not per window. τ is set by the first minimum of Average Mutual Information (Fraser & Swinney 1986, Phys Rev A 33:1134); m is set by Cao's method (Cao 1997, Physica D 110:43–50), bounded to m ≤ 5. These global estimates are then applied uniformly across all windows and all realisations to ensure a consistent phase-space reconstruction throughout the record.

**Recurrence threshold** is fixed by setting recurrence rate RR = 0.05 (Kraemer et al. 2018, Chaos 28:085720), making DET and TRANS directly comparable across windows and realisations regardless of signal amplitude variation.

**Theiler window** = max(1, τ(m − 1)), excluding autocorrelated neighbours along the line of identity (Theiler 1986; Marwan 2010).

Ensemble median and 5th/95th percentile envelopes across the 1,000 realisations propagate the full age-model uncertainty into both dynamical measures.

**Key configuration options** (top of script):

| Parameter | Default | Notes |
|-----------|---------|-------|
| `WINDOW_SIZE` | 100 | Sliding window width (points) |
| `WINDOW_STEP` | 5 | Step between windows (points) |
| `TARGET_RR` | 0.05 | Fixed recurrence rate |
| `MIN_LINE` | 2 | Minimum diagonal line length for DET |
| `MAX_M` | 5 | Upper bound on embedding dimension (Cao's method) |
| `SMOOTH_SIGMA` | 2.0 | Gaussian smoothing σ for display curves (points) |
| `DET_THRESHOLD` | 0.4 | Reference line and threshold for sustained-marker detection |
| `MIN_SUSTAIN` | 80 | Minimum consecutive windows below threshold to draw a marker |
| `N_REALISATIONS` | `None` (all 1,000) | Set e.g. `200` for a quick test run |
| `FORCE_TAU` / `FORCE_M` | `None` | Override AMI/Cao estimates with fixed integers |

A full 1,000-realisation run takes approximately 30–90 minutes depending on hardware.

---

## Data Availability

Raw proxy data (Co/Ca, Ni/Ca, δ¹⁸O depth series) and U/Th chronology for stalagmite HS4 (Heshang Cave, China) are included in `data/Drip_rate.xlsx`.

Full archive (data + outputs): [Zenodo DOI 10.5281/zenodo.16392750](https://doi.org/10.5281/zenodo.16392750)

---

## Funding and Acknowledgements

This research was supported by:

- **EU Horizon 2020** — Marie Skłodowska-Curie Actions, grant no. 691037 ([QUEST: QUantitative paleoEnvironments from SpeleoThems](https://quest.pik-potsdam.de/))
- **Te Apārangi — Royal Society of New Zealand**, grants RIS-UOW1501 and Rutherford Discovery Fellowship RDF-UOW1601
- **Ministry for Business, Innovation and Employment (MBIE)**, grant UOWX2102

Splash photograph © Garry Smith, from Hartland & Zitoun (2018), *Geochemical Perspectives Letters*.

---

## Authors and Affiliations

Adam Hartland¹², Bedartha Goswami³, Jungho Park¹, Sebastian Höpker², Dorisel Torres Rojas², Jin Liao⁴, Bethany R.S. Fox⁵, Norbert Marwan³, Sebastian F.M. Breitenbach⁶, Chaoyong Hu⁴

¹ Lincoln Agritech Ltd, Ruakura, Hamilton, Aotearoa New Zealand  
² Te Aka Mātuatua, School of Science, University of Waikato, Hamilton, Aotearoa New Zealand  
³ Potsdam Institute for Climate Impact Research, Potsdam, Germany  
⁴ China University of Geosciences, Wuhan, China  
⁵ School of Applied Sciences, University of Huddersfield, UK  
⁶ Department of Geography and Environmental Sciences, Northumbria University, UK  

Correspondence: [adam.hartland@lincolnagritech.co.nz](mailto:adam.hartland@lincolnagritech.co.nz)

---

## License

MIT License — see [LICENSE](LICENSE) for details.


