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

## 2. Install Python dependencies

Requires **Python 3.9+**. We recommend using a virtual environment or conda.

```bash
# Option A: pip (virtual environment)
python -m venv venv
source venv/bin/activate        # Linux/Mac
venv\Scripts\activate           # Windows

pip install flask numpy pandas scipy

# Option B: conda
conda create -n paleodriprates python=3.10 flask numpy pandas scipy
conda activate paleodriprates
```

### BayProX

The app depends on the `bayprox` package for Bayesian age–depth modelling. Install it from its own repository:

```bash
pip install bayprox
```

If `bayprox` is not pip-installable, clone it alongside `paleodriprates` and ensure it is on the Python path (the app adds `..` to `sys.path` at runtime).

---

## 3. Repository structure

```
paleodriprates/
├── dr_app/
│   ├── app.py                    ← Flask web application (main file)
│   ├── model.py                  ← Forward model: h(V) trace element kinetics
│   ├── params.py                 ← Physical constants (VMAX, VMIN, VRES, etc.)
│   ├── driprates_stochastic.py   ← Drip rate PDF computation with stochastic priors
│   ├── model_stochastic.py       ← Stochastic model wrapper
│   └── concentration_prior.py    ← Log-normal priors for [TE]aq and [Ca]aq
├── bayprox/                      ← BayProX age-depth engine (if not pip-installed)
│   ├── data.py
│   ├── agedepth.py
│   └── proxyrecord.py
└── README.md
```

---

## 4. Run the app

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

## 5. Prepare your data

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

## 6. Workflow

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

## 7. Troubleshooting

| Problem | Fix |
|---------|-----|
| Blank plots after run | Hard refresh: `Ctrl+Shift+R` (Windows) or `Cmd+Shift+R` (Mac). Check footer for build version. |
| `KeyError: 'Zn'` or similar | Set a numeric Kp value — theoretical Kp only supports Cu, Ni, Co. |
| Ages cluster at zero on plot | Hard refresh — this was fixed in v24+. |
| Heatmap appears featureless | Fixed in v39+ with per-column normalisation and V-range auto-crop. |
| Model stuck at a percentage | Check the terminal — look for Python tracebacks. Common cause: Kp=−1 for unsupported element. |
| Depth mismatch warning | Verify the depth unit dropdowns on both the Depth/Age and TE panels match your data. |

---

## 8. Citation

If you use this tool, please cite:

> Hartland, A., Goswami, B., Park, J., Höpker, S., Torres Rojas, D., Liao, J., Fox, B.R.S., Marwan, N., Breitenbach, S.F.M., Hu, C. *Quantitative Holocene precipitation reconstruction from stalagmite trace metal kinetics reveals East Asian monsoon drivers.* Submitted to Nature Geoscience.

Code & data: [github.com/waikatosci/paleodriprates](https://github.com/waikatosci/paleodriprates) · [doi:10.5281/zenodo.16392750](https://doi.org/10.5281/zenodo.16392750)

---

## 9. Contact

Adam Hartland — [adam.hartland@lincolnagritech.co.nz](mailto:adam.hartland@lincolnagritech.co.nz)

Lincoln Agritech Ltd / University of Waikato, New Zealand
