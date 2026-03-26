# PaleoDripRates — Quick Start Guide

**Drip Rate Estimator** is a browser-based tool for reconstructing past cave drip rates (and precipitation) from speleothem trace metal concentrations. It wraps the BayProX Bayesian age–depth modelling framework with a kinetic inversion engine based on organic–metal complex (OMC) dissociation.

---

## 1. Get the code

```bash
git clone https://github.com/waikatosci/paleodriprates.git
cd paleodriprates
```

Or download and extract the [Zenodo archive](https://doi.org/10.5281/zenodo.16392750).

---

## 2. Install Python

Requires **Python 3.10+** (3.12 recommended). If you already have Python installed, skip to step 3.

**Windows:**
1. Download from https://www.python.org/downloads/
2. Run the installer — **tick "Add python.exe to PATH"** on the first screen before clicking Install
3. Restart any open command prompts
4. Verify: `python --version`

**Mac/Linux:** Python 3 is typically pre-installed. Verify with `python3 --version`. If not present, install via Homebrew (`brew install python`) or your package manager.

---

## 3. Install Python dependencies

```bash
pip install flask numpy pandas scipy
```

That's the only pip command needed. Everything else is included in the repository.

Optionally, use a virtual environment to keep dependencies isolated:

```bash
# Option A: pip (virtual environment)
python -m venv venv
source venv/bin/activate        # Linux/Mac
venv\Scripts\activate           # Windows
pip install flask numpy pandas scipy

# Option B: conda
conda create -n paleodriprates python=3.12 flask numpy pandas scipy
conda activate paleodriprates
```

### BayProX

The app depends on the `bayprox` package for Bayesian age–depth modelling. Install it from its own repository:

```bash
pip install bayprox
```

If `bayprox` is not pip-installable, clone it alongside `paleodriprates` and ensure it is on the Python path (the app adds `..` to `sys.path` at runtime).

---

## 4. Repository structure

```
paleodriprates/
├── dr_app/
│   ├── app.py                    ← Flask web application (main file)
│   ├── model.py                  ← Forward model: h(V) trace element kinetics
│   ├── params.py                 ← Physical constants (VMAX, VMIN, VRES, etc.)
│   ├── driprates_stochastic.py   ← Drip rate PDF computation with stochastic priors
│   ├── model_stochastic.py       ← Stochastic model wrapper
│   ├── concentration_prior.py    ← Log-normal priors for [TE]aq and [Ca]aq
│   └── HS4_example_inputs/       ← Example data (Heshang Cave stalagmite HS4)
│       ├── HS4_depth_age.csv     ← U-Th dating table
│       └── HS4_trace_elements.csv← Co and Ni LA-ICP-MS profiles
├── bayprox/                      ← BayProX age-depth engine (if not pip-installed)
│   ├── data.py
│   ├── agedepth.py
│   └── proxyrecord.py
├── QUICKSTART.md                 ← This file
└── README.md
```

---

## 5. Run the app

```bash
cd paleodriprates/dr_app
python app.py
```

You should see:

```
 ===================================
  Drip Rate Estimator - Starting...
 ===================================

  Open your browser at:  http://localhost:5000
```

Open **http://localhost:5000** in any modern browser (Chrome, Firefox, Edge).

---

## 6. Quick test with HS4 example data

The repository includes example input files from **stalagmite HS4, Heshang Cave, China** — the validation dataset from the manuscript. This lets you verify your installation produces correct output in under 15 minutes.

### Step 1 — Upload data

1. Navigate to **Data Inputs**
2. Upload `dr_app/HS4_example_inputs/HS4_depth_age.csv` to the **Depth / Age** dropzone
3. Upload `dr_app/HS4_example_inputs/HS4_trace_elements.csv` to the **Trace Elements** dropzone
4. The app will auto-detect columns and elements — verify the age–depth plot looks reasonable (monotonically increasing, 0–260 cm depth, 0–10,000 yrs BP)

### Step 2 — Accept defaults and run

All model parameters are pre-configured with sensible defaults for HS4:
- Cave temperature, Ca concentration (66.46 ppm), and element-specific Kp/Kd values auto-populate
- Aqueous TE concentrations from Heshang Cave monitoring: Co = 0.47 ppb, Ni = 4.34 ppb
- Concentrations are fixed (no stochastic priors) by default — this produces the tightest, most physically interpretable reconstruction
- The element selector detects Co and Ni from the column headers

Simply navigate to **Run Model** and click **Run Model**. Typical runtime is 10–15 minutes for two trace elements over ~10,000 years.

### Step 3 — Check outputs

When complete, click **Results** to explore the four output tabs:

**📈 Time Series** — reconstructed drip rate showing the Holocene decline from ~20 drips/min in the early Holocene, with the pronounced 5.2 ka aridity event clearly resolved as a sharp dip, and the modern (post-1820) increase.

**🌈 PDF Heatmap** — full probability density field showing how uncertainty varies through time. Select the **Reds** colourmap for clearest contrast. Toggle **Percentiles** to overlay the median and IQR bands.

**💧 Smart & Friedrich** — hydrological classification placing HS4 in the **Buffered overflow** quadrant (mean ~20 drips/min, CV ~0.25), consistent with the site's known hydrology.

**📐 Age Model** — BayProX posterior age–depth relationship with dated points, showing the smooth monotonic age model through 260 cm of growth.

All outputs are downloadable from the **Downloads** tab as CSV/JSON files.

---

## 7. Prepare your own data

You need **two CSV files**:

### Depth–Age CSV

| Column | Description |
|--------|-------------|
| `Depth` | Sample depth (cm or mm) |
| `Age(yr BP)` | Calendar age (years before present) |
| `error` | 2σ age uncertainty |

### Trace Element CSV

| Column | Description |
|--------|-------------|
| `Depth` or `Depth(mm)` | Measurement depth — unit auto-detected from header name |
| One or more TE columns | e.g. `Cu_ppm`, `Zn_ppm`, `Co_ppm` — element auto-detected |

**Depth units:** The app auto-detects `mm` or `cm` from column names (e.g. `Depth(mm)` → mm). If both files use different units, the backend automatically converts TE depths to match the depth–age file. You can also override with the unit dropdown.

---

## 8. Full workflow

### Step 1 — Data Inputs

1. Upload the **Depth / Age** CSV and verify column assignments
2. Upload the **Trace Elements** CSV
3. Select the depth column and one or more proxy columns
4. Review the age–depth plot and TE preview chart
5. If your record has a hiatus, use the **Growth Rate & Hiatus Detection** panel to auto-detect or manually define exclusion zones

### Step 2 — Model Parameters

1. Set **cave temperature** and **monitored drip rate**
2. Set **Ca aqueous concentration** (manual or from monitoring CSV)
3. For each TE, review:
   - **Element** — auto-detected from column name
   - **Kp** — partition coefficient (auto-filled from literature for supported elements; Wang & Xu theoretical for Cu/Ni/Co)
   - **ln(Kd)** — dissociation rate constant (auto-calculated from drip rate when in drip-rate mode)
   - **Calibration window** — % of recent deposition used to anchor Kd
4. Check the **h(V) diagnostic chart** — your observed data should fall within the model curve

> **Note on concentration uncertainty:** The default run uses fixed mean concentrations for Ca and TE aqueous phases. The Std dev fields beside each concentration input are intentionally left blank. Entering standard deviations activates stochastic marginalisation over log-normal concentration priors, which can substantially inflate uncertainty envelopes and shift the median reconstruction (see Supplementary Information). The fixed-concentration approach is recommended for standard use; stochastic priors should only be explored as a sensitivity test.

### Step 3 — Run Model

1. Review the **pre-flight summary** and sanity checks
2. Click **Run Model** — typical runtime is 5–15 minutes for two TEs over ~12,000 years
3. Progress is shown in the status bar

### Step 4 — Results

- **📈 Time Series** — median drip rate with IQR and 5–95% envelopes
- **💧 Smart & Friedrich** — hydrological classification (seepage / fracture / overflow)
- **🌈 PDF Heatmap** — full probability density with percentile overlay
- **📐 Age Model** — BayProX posterior age–depth relationship

Use the **Drip rate / τ toggle** to switch between drip rate (min⁻¹) and thin-film residence time τ (s) for non-stalagmite applications.

### Downloads

All outputs are available as CSV/JSON files from the **Downloads** tab, including:
- `drip_rate_summary.csv` — percentiles at each timestep
- `drip_rate_realisations.csv` — full ensemble for RQA
- `age_model.csv` — BayProX age model
- `input_summary.csv` — complete parameter record for reproducibility

---

## 9. Troubleshooting

| Problem | Fix |
|---------|-----|
| Blank plots after run | Hard refresh: `Ctrl+Shift+R` (Windows) or `Cmd+Shift+R` (Mac). Check footer for build version. |
| `KeyError: 'Zn'` or similar | Set a numeric Kp value — theoretical Kp only supports Cu, Ni, Co. |
| Ages cluster at zero on plot | Hard refresh — this was fixed in v24+. |
| Heatmap appears featureless | Fixed in v39+ with per-column normalisation and V-range auto-crop. |
| Model stuck at a percentage | Check the terminal — look for Python tracebacks. Common cause: Kp=−1 for unsupported element. |
| Depth mismatch warning | Verify the depth unit dropdowns on both the Depth/Age and TE panels match your data. |

---

## 10. Citation

If you use this tool, please cite:

> Hartland, A., Goswami, B., Park, J., Höpker, S., Torres Rojas, D., Liao, J., Fox, B.R.S., Marwan, N., Breitenbach, S.F.M., Hu, C. *Quantitative Holocene precipitation reconstruction from stalagmite trace metal kinetics reveals East Asian monsoon drivers.* Submitted to Nature Geoscience.

Code & data: [github.com/waikatosci/paleodriprates](https://github.com/waikatosci/paleodriprates) · [doi:10.5281/zenodo.16392750](https://doi.org/10.5281/zenodo.16392750)

---

## 11. Contact

Adam Hartland — [adam.hartland@lincolnagritech.co.nz](mailto:adam.hartland@lincolnagritech.co.nz)

Lincoln Agritech Ltd / University of Waikato, New Zealand
