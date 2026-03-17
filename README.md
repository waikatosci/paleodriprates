# PaleodripRates: Kinetic Proxy for Stalagmite Drip Rate and Precipitation Reconstruction

This repository contains code, data, and utilities for reconstructing cave drip rates and Holocene precipitation from stalagmite trace metals using a kinetic proxy model based on organic-metal complex dissociation. The methods are detailed in:

> Hartland, A., Goswami, B., Höpker, S.N., Park, J., Torres Rojas, D., Liao, J., Fox, B.R.S., Marwan, N., Breitenbach, S.F.M., & Hu, C. (2025). Quantitative Holocene precipitation reconstruction from stalagmite trace metal kinetics reveals East Asian Monsoon drivers. *Nature Geoscience*. DOI: [insert upon publication]

The proxy exploits the distributed dissociation kinetics of transition metals (e.g., Co, Ni) bound to organic ligands in dripwater, calibrated against modern data from Heshang Cave, China. Drip rates are inverted probabilistically via Bayesian methods, then chained to precipitation/temperature (P/T) regressions and temperature proxies for quantitative precipitation estimates (mm yr⁻¹) with propagated uncertainties.

---

## Contents

- [Dependencies](#dependencies)
- [Installation](#installation)
- [Web Application (app.py)](#web-application-apppy)
  - [Launching the app](#launching-the-app)
  - [Data Input](#data-input)
  - [Model Parameters](#model-parameters)
  - [Analysis Modes](#analysis-modes)
  - [Output Options](#output-options)
  - [Running the model](#running-the-model)
  - [Output Explorer](#output-explorer)
- [Script-based Usage](#script-based-usage)
- [Repository Structure](#repository-structure)
- [Data Availability](#data-availability)
- [Citation](#citation)
- [License](#license)
- [Acknowledgments](#acknowledgments)

---

## Dependencies

- Python 3.8+
- numpy, pandas, scipy (integrate, optimize, interpolate, stats)
- flask
- matplotlib
- openpyxl, Pillow (PIL), progressbar
- bayprox (custom Bayesian library included under `lib/`)

---

## Installation

```bash
git clone https://github.com/waikatosci/paleodriprates.git
cd paleodriprates
pip install -r requirements.txt
```

Or install essential components directly:

```bash
pip install numpy pandas matplotlib scipy flask openpyxl pillow progressbar
```

`bayprox` is a custom module included under `lib/`. Add it to your `PYTHONPATH` if importing outside the repo root.

---

## Web Application (app.py)

`dr_app/app.py` is a self-contained Flask web application providing a browser-based interface to the full reconstruction pipeline. It replaces the legacy Excel-based workflow.

### Launching the app

**macOS / Linux:**
```bash
bash launch_mac_linux.sh
```

**Windows:**
Double-click `launch_windows.bat` (requires Anaconda or a Python environment with Flask).

Then open your browser at **http://localhost:5000**. A hard-refresh (`Ctrl+Shift+R`) is recommended after updating `app.py` to clear any cached pages.

---

### Data Input

The Data Input panel accepts CSV files with header rows. All depth values should be in **cm**.

#### Depth / Age Model (required)

Upload a CSV containing depth, calibrated age (yrs BP), and age error (2σ) columns. After upload, column selectors appear and an **age–depth plot** is rendered automatically showing:

- Dated control points as scatter markers with tooltips (depth → age).
- A **fit curve** passing through all dated points, with the portion within the data range drawn as a solid green line and extrapolated portions as dashed grey lines.
- A shaded band indicating the currently selected calibration age window.

A **Fit** selector in the plot header switches between two models:

| Mode | Behaviour | When to use |
|------|-----------|-------------|
| **Monotone spline (PCHIP)** *(default)* | Piecewise Cubic Hermite Interpolating Polynomial (Fritsch-Carlson 1980). Passes exactly through every dated point. Enforces monotonicity — age always increases with depth, no oscillations. Handles kinks in the accumulation rate (e.g. slow-growth intervals) correctly. | Most records, especially those with variable or interrupted accumulation rates. |
| **Linear regression** | Ordinary least-squares fit through all points. Does not pass through individual points. | Simple records with uniform accumulation; useful as a sanity check. |

Extrapolation beyond the data range uses the end tangents of the spline (linear extension from the first and last knot slopes), avoiding unbounded cubic behaviour outside the dated interval. The curve is sampled at 300 evenly spaced points for a smooth visual.

The **Calibration age min/max** fields in Site Settings are auto-populated from the extrapolated fit endpoints on upload and update when the column selection or fit type changes. They can be overridden at any time by editing the extrap handles below the plot or the fields directly; both stay in sync.

#### Trace Elements (required)

A single CSV upload accepts one or more trace element proxy columns alongside a depth column. After upload:

- A **depth column** selector appears.
- One or more **proxy rows** are shown, each with a column selector and a concentration unit selector (ppb, ppm, µg/g, mg/kg). The unit selection triggers automatic conversion to ppb in the backend before the model runs.
- Element identity is **auto-detected** from the column name (e.g. a column named `Co_ppb`, `Cobalt`, or `co` will pre-select Co in the corresponding parameter card). Detection fires on upload, on column change, and when rows are added.
- The **+ Add trace element** button appends a new row. Rows can be removed with the × button (minimum one row required).
- Multiple elements are combined via geometric mean of their individual drip rate PDFs.

#### Data Preprocessing (optional but recommended for high-resolution data)

After upload a **Data Preprocessing** panel appears below the proxy row selectors. LA-ICP-MS and other high-resolution datasets (e.g. 1500+ rows at 10–100 µm spacing) can substantially increase BayProX runtime because the unified depth grid scales with the number of depth points. The preprocessing panel provides two tools applied server-side before the model runs, so the saved `trace_elem1.csv` reflects the cleaned data:

**Sigma-clip outlier removal** — each numeric proxy column is inspected independently using a centred rolling window. For each point, the local mean and standard deviation are computed from the surrounding *N* neighbours (default window = 50 points). Points deviating more than *σ* standard deviations from their local mean are flagged as outliers and excluded from the block average. A window of 0 falls back to global sigma-clipping (single mean and SD across the whole record), which is simpler but can be too conservative for long records with a real secular trend — the rolling window approach avoids this by testing each point against its local neighbourhood rather than the record-wide statistics. This is the same algorithm used by the backend `detect_outliers` / `clean()` pipeline, so the preview faithfully represents what the model will see.

A mini scatter plot shows kept points (green) and flagged outliers (red) in real time as the threshold or window size is adjusted, without committing any changes until **Apply & Save** is clicked.

**Block-average downsampling** — the surviving rows are grouped into blocks of size `⌊N_rows / target_N⌋` and each block is replaced by its column-wise mean. This preserves the signal shape while reducing computational cost. The target can be set as an absolute row count. The before/after badge above the chart updates as the target is adjusted.

Clicking **Apply & Save** sends both parameters to the `/preprocess` endpoint, which overwrites `trace_elem1.csv` and returns the new row count. The preprocessing status line confirms the result and the runtime estimate (see below) recalculates automatically.

> **Guidance for LA-ICP-MS data:** typical BayProX runs complete in 20–40 minutes with ~300–500 depth points. If your data are already averaged to 100 µm spacing, block-averaging to ~300 points is generally sufficient to retain all palaeoclimate-relevant variability while keeping runtimes tractable. The default sigma-clip window of 50 points and threshold of 3σ is conservative; for data with known instrumental spikes a window of 20–30 points and threshold of 2.5σ may be more appropriate. Set window size to 0 to use global sigma-clipping if the record lacks a long-term trend.

#### Isotope Data (optional)

Upload a CSV with depth and isotope (e.g. δ¹⁸O in ‰) columns. This dataset is used to build the unified depth grid and proxy record but does not directly enter the drip rate PDF calculation. If omitted, the unified depth grid is built from the trace element depths alone.

---

### Model Parameters

#### Analysis Mode

A toggle at the top of the parameters panel selects between two operational modes. See [Analysis Modes](#analysis-modes) below.

#### Element guidance

Elements exhibiting OMC dissociation kinetics suitable for drip rate reconstruction are primarily d-block transition metals. Redox-active elements (Mn, Fe) are not suitable. Elements beyond Ni, Co, and Cu are not yet characterised in cave systems. Kp = −1 instructs the model to calculate the partition coefficient theoretically from cave temperature using the Wang & Xu (2001) lattice strain model; positive values override this.

#### Per-element parameter cards

One card is shown for each configured trace element row. Parameters include:

| Parameter | Description |
|-----------|-------------|
| Element | Selects element and auto-fills molecular weight and default Kp |
| Kp | Partition coefficient; −1 = theoretical (Wang & Xu 2001) |
| Mean ln(Kd) | Log-mean of the Kd distribution |
| Std dev ln(Kd) | Log-standard deviation of Kd |
| XF — Fast fraction | Labile fraction dissociating too rapidly to fractionate with drip rate |
| XI — Inert fraction | Fraction not participating in OMC dissociation |
| XL — Labile fraction | Auto-calculated as 1 − XI − XF |
| Aqueous concentration | Element concentration in drip water (unit selectable) |

For Ni the recommended defaults are XI = 0.10, XF = 0.01, XL = 0.89 (Kd mean/sd from Lindeman et al., *GCA* 317, 2022).

Cards for TE3 and beyond are generated dynamically when rows are added in Data Input and are hidden when not needed.

#### Cave Conditions

Cave temperature (°C) and Ca concentration in drip water (unit selectable). These are used in full quantification mode only; in semi-quantitative mode they are hidden.

---

### Analysis Modes

#### Full Quantification (default)

Requires aqueous element concentrations, Ca concentration, and calibrated Kd / Kp values. Outputs absolute drip rate reconstructions in **drips min⁻¹** suitable for hydrological quantification and Smart & Friedrich classification.

#### Semi-Quantitative

Aqueous chemistry inputs (element and Ca concentrations) are hidden. Element fractions and Kd values still control the *shape* of the drip rate response and should be set from the literature. Output is normalised to a reference value and expressed as **% of reference drip rate**.

Three normalisation options are available (evaluated in priority order):

1. **Anchor drip rate** — if a single modern monitoring measurement is available (drips min⁻¹), enter it here. This converts the normalised output back to absolute units without requiring full water chemistry.
2. **Reference period** — specify an age window (yrs BP). The mean of the pc50 time series within that window becomes 100%. Useful when modern calibration data exist for a specific interval.
3. **Record maximum** (default fallback) — the 95th percentile maximum across the record is set to 100%.

In semi-quantitative mode the realisations ensemble (if generated) is normalised independently per realisation, preserving the full distribution of uncertainty through the normalisation step — so RQA analyses on the realisations remain statistically consistent with the summary percentiles.

The Smart & Friedrich classification tab is disabled in semi-quantitative mode as the absolute mean discharge is not defined.

---

### Output Options

#### Realisations (for RQA)

A toggle controls whether the full ensemble of drip rate realisations is generated and saved. When enabled, `n_realisations` independent draws are sampled from the joint posterior PDF at each time step. The output CSV has shape *(N_timesteps + 1) × N_realisations* and is suitable for direct input to recurrence quantification analysis (RQA) tools. Random seed is configurable for reproducibility.

When disabled, only the summary percentile CSV is written, which substantially reduces computation time for exploratory runs.

#### Proxy Record caching

The BayProX proxy record computation is the slowest step (progress and ETA are shown on the Run page). If you have already run the pipeline once and are only changing model parameters (Kd, Kp, fractions, concentrations), check **Use cached proxy record** to skip straight to the drip rate step using the saved `ProxyRecord.pkl`.

---

### Running the model

The Run page shows a **Configuration Summary** auto-populated from current inputs, reflecting the actual number of configured trace elements, their column names and units, and the selected analysis mode.

#### Runtime estimate

Before starting, a **Runtime Estimate** card shows three statistics derived from the current configuration:

| Stat | Source |
|------|--------|
| Depth points | Row count of the (preprocessed) TE CSV |
| Age range | `calage_max − calage_min` from Site Settings |
| Est. runtime | `depth_points × age_range × n_TEs × empirical_constant + overhead` |

The empirical constant (`2.8 × 10⁻⁵` s per depth-point per age-step) was calibrated against typical laptop hardware running BayProX. Runtimes on faster machines will be lower; on older hardware they may be higher. The estimate turns amber when the projected runtime exceeds 20 minutes, and a contextual tip suggests either using the preprocessing panel to downsample or enabling **Use cached proxy record** if BayProX has already run for this dataset.

The estimate recalculates automatically whenever the TE CSV is uploaded, preprocessing is applied, the calage range changes, or the Run panel is opened.

#### Progress and ETA

A progress bar and percentage track pipeline stages. The **ETA label** cycles through three states:

- *Starting…* for the first 2 seconds.
- *Elapsed: X* until 5% progress is reached (insufficient data for a stable estimate).
- *ETA ~X* from 5% onwards based on linear extrapolation of elapsed time vs progress.
- *Completed in X* when the run finishes.

The log box streams backend messages in real time. Errors are highlighted in red. On completion the **View Results →** button becomes available.

---

### Output Explorer

The results panel is a tabbed explorer with two views.

#### Time Series

Reconstructed drip rate (or relative drip rate in semi-quantitative mode) over time with a full uncertainty envelope:

- 5–95th percentile shaded outer band.
- 25–75th percentile (IQR) inner band.
- Median line.
- **Kd sensitivity ribbon** — two amber dashed lines showing the median reconstruction at Kd_mn − Kd_sd and Kd_mn + Kd_sd (near-deterministic passes with Kd_sd ≈ 0). The gap between the amber ribbon and the green envelope reflects how much reconstruction uncertainty is attributable to partition coefficient literature uncertainty vs the proxy record itself.

The y-axis label switches between "Drip rate (min⁻¹)" (full mode) and "Relative drip rate (% of reference)" (semi-quantitative mode).

#### Smart & Friedrich (1987) Classification

Available in full quantification mode only. An x-y scatter plot of mean discharge (y, drips min⁻¹) vs coefficient of variation (x, CV = σ/μ of the median time series) with four coloured classification zones:

| Zone | CV | Mean discharge | Interpretation |
|------|----|----------------|----------------|
| Seepage / percolation | Low (< 0.5) | Low (< 5 min⁻¹) | Slow matrix flow; highly attenuated signal |
| Fracture / conduit | High (≥ 0.5) | Low | Episodic fracture recharge; strong seasonal signal |
| Buffered overflow | Low | High (≥ 5 min⁻¹) | Perched water table; sustained but modulated |
| Flood / conduit overflow | High | High | Rapid conduit response to recharge events |

The reconstruction point is plotted with IQR error bars on both axes (mean: temporal Q25–Q75 of the pc50 series; CV: propagated from the pc05/pc95 spread). A classification label names the zone.

---

### Output files

| File | Contents |
|------|----------|
| `drip_rate_summary.csv` | Percentile columns `age, pc05, pc10, pc25, pc50, pc75, pc90, pc95` at each time step |
| `drip_rate_realisations.csv` | Full ensemble; shape *(N_timesteps + 1) × N_realisations*; suitable for RQA |
| `ProxyRecord.pkl` | Cached BayProX proxy distributions; reused if *Use cached proxy record* is checked |
| `uploads/trace_elem1.csv` | Overwritten in-place when preprocessing is applied; contains the block-averaged, sigma-clipped data that actually enters the pipeline |

---

## Script-based Usage

The original script-based pipeline remains available for batch processing and headless execution.

**1 — Kinetic model and dissociation kinetics**

Run `model.py` for numerical integration of dissociated fractions using a log-normal distribution of rate constants. Example: compute the labile fraction for a given residence time τ.

**2 — Bayesian drip rate inversion**

Use `drip_rate_parallel.py` (recommended) or `drip_rate_serial.py` (slower). Input: `Drip_rate.xlsx` containing depth-age data and proxy records (Co/Ni concentrations) processed via BayProX for age-depth uncertainties. Output: drip rate PDFs (medians, percentiles) written to the `PlotDripRate` and `PlotDripIsotope` sheets.

**3 — Precipitation reconstruction**

Run `P_quantification_Holocene.ipynb` for chained P/T regressions and Monte Carlo propagation. Inputs: drip rate percentiles and temperature estimates (`Precip_from_drip_rates.xlsx`). Outputs: `p_reconstruction.csv` (annual precipitation posteriors with medians and 25–75th percentiles) and `p_plot.png`. See `precip_recon/readme.txt` for supplementary details.

**4 — Utilities**

`utils.py` provides progress bars for long computations. Example data (`Drip_rate.xlsx`, `ProxyRecordPlot.xlsx`) are included for calibration and sensitivity tests.

---

## Repository Structure

```
paleodriprates/
├── dr_app/
│   ├── app.py                  # Flask web application (all pipeline stages)
│   │                           #   Routes: / · /upload · /preprocess · /run
│   │                           #           /status · /chart_data · /download
│   ├── drip_rate_util.py       # Bayesian drip rate inversion utilities
│   ├── params.py               # Outlier detection parameters
│   ├── launch_mac_linux.sh     # macOS/Linux launcher
│   ├── launch_windows.bat      # Windows launcher
│   └── uploads/                # Runtime CSV upload directory (auto-created)
│   └── outputs/                # Runtime output directory (auto-created)
├── model.py                    # Core dissociation kinetics and expectation calculations
├── drip_rate_parallel.py       # Parallel Bayesian drip rate inversion (script)
├── drip_rate_serial.py         # Serial version
├── P_quantification_Holocene.ipynb  # Precipitation reconstruction notebook
├── utils.py                    # Helper functions (progress bars, etc.)
├── lib/
│   └── bayprox/                # Custom Bayesian proxy library
├── data/
│   ├── Drip_rate.xlsx          # Depth-age and proxy input data
│   ├── Precip_from_drip_rates.xlsx
│   └── ProxyRecordPlot.xlsx    # Calibration and sensitivity data
├── precip_recon/
│   └── readme.txt
├── figures/                    # Generated plots
├── requirements.txt
└── LICENSE
```

---

## Data Availability

Raw proxy data and monitoring records are included in the Excel files in `data/`. Full archive (including all calibration datasets): **Zenodo DOI: 10.5281/zenodo.16392750**.

---

## Citation

**Manuscript:**
> Hartland, A., Goswami, B., Höpker, S.N., Park, J., Torres Rojas, D., Liao, J., Fox, B.R.S., Marwan, N., Breitenbach, S.F.M., & Hu, C. (2025). Quantitative Holocene precipitation reconstruction from stalagmite trace metal kinetics reveals East Asian Monsoon drivers. *Nature Geoscience*. DOI: [insert upon publication]

**Repository:**
> Hartland, A. et al. (2025). PaleodripRates: Code for stalagmite drip rate and precipitation reconstruction. Zenodo. https://doi.org/10.5281/zenodo.16392750

---

## License

MIT License — see `LICENSE` for details.

---

## Acknowledgments

Funded by EU Horizon 2020 Marie Skłodowska-Curie Actions (no. 691037, QUEST), Te Apārangi Royal Society of New Zealand (RIS-UOW1501), Ministry for Business, Innovation and Employment (UOWX2102), and a Rutherford Discovery Fellowship (RDF-UOW1601).

For questions or contributions, open an issue on GitHub or contact the corresponding author: adam.hartland@lincolnagritech.co.nz
