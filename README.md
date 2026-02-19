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


Usage
1) Kinetic Model and Dissociation Kinetics:
Run model.py for numerical integration of dissociated fractions (e.g., log-normal distribution of rate constants).
Example: Compute labile fraction for given residence time (τ).

2) Bayesian Inversion for Estimating Drip Rates:
Use drip_rate_parallel.py (recommended for faster performance) or drip_rate_serial.py (slower) to estimate drip rates from trace metal posteriors.
Input: Drip_rate.xlsx: Depth-age data, Proxy records (Co/Ni concentrations) processed via BayProx for age-depth uncertainties, and oxygen isotope data.
Output: Drip rate PDFs (medians, percentiles) saved in the PlotDripRate and PlotDripIsotope sheets of drip_rate.xlsx.

3) Precipitation Reconstruction:
Run P_quantification_Holocene.ipynb (Jupyter notebook) for chained regressions and Monte Carlo propagation.
Inputs: Drip rate percentiles, temperature estimates (Precip_from_drip_rates.xlsx).
Outputs: Annual precipitation posteriors (p_reconstruction.csv) with medians and 25–75th percentiles.
Visualization: Generates plots like p_plot.png for Holocene P trends.
Supplementary: See readme.txt file in precip_recon directory

5) Utilities:
utils.py: Progress bars for long computations.
Example data: Drip_rate.xlsx, ProxyRecordPlot.xlsx for calibration and sensitivity tests.
For executable versions integrated with Excel, contact the corresponding author (adam.hartland@lincolnagritech.co.nz or waikatoscientific@gmail.com). All scripts are self-contained; test with provided snippets in the manuscript supplement.

Repository Structure
model.py: Core dissociation kinetics and expectation calculations.
drip_rate_util.py: Bayesian inversion for drip rates.
P_quantification_Holocene.ipynb: Notebook for precipitation reconstruction via Monte Carlo.
utils.py: Helper functions (e.g., progress bars).
lib/: Custom libraries (e.g., bayprox).
data/: Excel files for inputs/outputs (e.g., drip_rate_percentiles.xlsx, Precip_from_drip_rates.xlsx).
figures/: Generated plots (e.g., sensitivity analyses, reconstructions).
requirements.txt: List of dependencies.
LICENSE: MIT License (or as specified).

Data Availability
Raw proxy data and monitoring records are included in Excel files.
Full archive on Zenodo: DOI: 10.5281/zenodo.16392750.
Upon publication, an executable Excel-integrated version will be added.

Citation
If using this code or data, please cite:
Hartland, A., Goswami, B., Höpker, S.N., Park, J., Torres Rojas, D., Liao, J., Fox, B. R. S., Marwan, N., Breitenbach, S. F. M., & Hu, C. (2025). Quantitative Holocene precipitation reconstruction from stalagmite trace metal kinetics reveals East Asian monsoon drivers. Nature Geoscience. DOI: [insert DOI upon publication].

For the repository:
Hartland, A. et al. (2025). PaleodripRates: Code for stalagmite drip rate and precipitation reconstruction. Zenodo. https://doi.org/10.5281/zenodo.16392750

License
This project is licensed under the MIT License—see the LICENSE file for details.

Acknowledgments
Funded by EU Horizon 2020 Marie Skłodowska-Curie (no. 691037, QUEST), Te Apārangi Royal Society of New Zealand (RIS-UOW1501), Ministry for Business, Innovation and Employment (UOWX2102), and Rutherford Discovery Fellowship (RDF-UOW1601).

For questions or contributions, open an issue or contact the corresponding author.
