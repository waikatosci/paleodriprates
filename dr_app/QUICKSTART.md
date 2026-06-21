# Dr Paleo — Quick Start Guide

Dr Paleo (dr_paleo = drip rate) is a browser-based tool for reconstructing past cave drip rates and precipitation from speleothem trace metal concentrations. It wraps the BayProX Bayesian age-depth modelling framework with a kinetic inversion engine based on organic-metal complex (OMC) dissociation.

---

## 1. Get the code

```bash
git clone https://github.com/waikatosci/paleodriprates.git
cd paleodriprates
```

Or download and extract the [Zenodo archive](https://doi.org/10.5281/zenodo.16392750).

---

## 2. Install Python dependencies

Requires Python 3.9 or later. We recommend a virtual environment or conda.

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

`bayprox` is bundled with this repository, so no separate installation is needed. When you run the app from `dr_app/`, it adds the repo root to the Python path at runtime so the bundled package is found automatically.

---

## 3. Launch Dr Paleo

```bash
cd dr_app
python app.py
```

You should see:

```
  Dr Paleo — Drip Rate Estimator
  Open your browser at:  http://localhost:5000
```

Open `http://localhost:5000` in any modern browser (Chrome, Firefox, Edge).

- Windows shortcut: double-click `launch_windows.bat`
- macOS/Linux shortcut: run `./launch_mac_linux.sh`

> A terminal window stays open while Dr Paleo is running. Leave it open; close it when you are done to shut down the server.

---

## 4. Quick test with HS4 example data

The repository includes example input files from stalagmite HS4, Heshang Cave, China, the validation dataset from the manuscript. This lets you verify your installation produces correct output in under 15 minutes.

### Step 1 — Upload data

1. Navigate to Data Inputs.
2. Upload `dr_app/HS4_example_inputs/HS4_age_depth.csv` to the Depth / Age dropzone.
3. Upload `dr_app/HS4_example_inputs/HS4_TE.csv` to the Trace Elements dropzone.
4. Dr Paleo auto-detects the columns; verify the dropdowns map depth, age, and proxy columns correctly.

### Step 2 — Check parameters

1. Navigate to Model Parameters.
2. All defaults are pre-set for HS4; you should not need to change anything.
3. Verify that the element selectors show Co (TE1) and Ni (TE2).
4. The monitored drip rate, Kp, Kd, and kinetic fractions are all populated from the literature defaults.

### Step 3 — Run

1. Navigate to Run.
2. Click Run Model.
3. Dr Paleo shows a progress bar and live log. ETA appears after about 5% progress. Typical runtime with HS4 defaults is 5-15 minutes depending on your machine.

### Step 4 — Check outputs

1. Navigate to Results when Dr Paleo reports completion.
2. You should see:
   - Drip rate vs time: time series with percentile envelope and Kd sensitivity ribbon
   - Probability density heatmap: column-normalised PDF across the V grid
   - Smart & Friedrich classification: hydrological regime zones
   - Age model: BayProX posterior with dating points and error bars
3. Download the output CSVs from the Downloads panel.

### Expected outputs

| File | Description |
|------|-------------|
| `drip_rate_summary.csv` | Percentile summary (pc05, pc25, pc50, pc75, pc95) at each timestep |
| `drip_rate_realisations.csv` | Full MC ensemble (age × N realisations) for RQA |
| `age_model.csv` | Depth-age mapping with errors |
| `input_summary.csv` | All input parameters for reproducibility |

---

## 5. Using your own data

### Required files

| File | Format | Required columns |
|------|--------|-----------------|
| Depth / Age | CSV | depth (cm), age (years BP), age error |
| Trace Elements | CSV | depth (cm or mm, auto-detected), one or more proxy columns |

### Optional files

| File | Purpose |
|------|---------|
| Isotope | δ¹⁸O or δ¹³C for joint proxy reconstruction |
| Aqueous monitoring | TE and Ca time series for stochastic concentration priors |

### Workflow

1. Upload your CSVs in the Data Inputs panel.
2. Map columns using the dropdowns; Dr Paleo auto-detects common naming patterns.
3. Select elements: choose the trace element for each TE row. Dr Paleo populates Kp, Kd, kinetic fractions, and aqueous concentration defaults from the literature.
4. Set cave conditions: temperature (°C), Ca concentration, and monitored drip rate.
5. Choose analysis mode:
   - Full quantification: absolute drip rate (drips min⁻¹). Requires solution chemistry.
   - Semi-quantitative: relative drip rate (% of reference). No water chemistry needed.
6. Run. Dr Paleo shows progress and ETA.
7. Explore results: interactive charts and downloadable CSVs.

### Tips

- Preprocessing: use the built-in sigma-clip and block-average tools (Data Inputs panel) to reduce noisy high-resolution data before running.
- Caching: after the first run, tick "Use cached proxy record" to skip the BayProX step on subsequent runs and go straight to drip rate computation.
- Hiatus detection: specify hiatus depths in the Advanced section if your stalagmite has known growth hiatuses.
- Kp = −1: enter −1 for Kp to use the theoretical/empirical default. Dr Paleo resolves this via Wang & Xu (2001) for Cu/Ni/Co, Lindeman (2022) empirical values, or literature defaults for other elements.
- NaN handling: missing values in your CSV are handled automatically; they appear as gaps in the reconstruction, not errors.

---

## 6. Troubleshooting

| Symptom | Cause / Fix |
|---------|-------------|
| Browser shows old UI after update | Hard refresh: `Ctrl+Shift+R` (Windows/Linux) or `Cmd+Shift+R` (Mac) |
| Upload fails with a JSON parse error | Your CSV likely has missing values; the app handles these, so hard-refresh and re-upload |
| Model run reports an error | Check the terminal window for a Python traceback |
| ETA keeps climbing | Large datasets (>1000 rows) are slow; use the preprocessing tools to downsample |
| Heatmap appears featureless | Check the V-range and per-column normalisation settings |
| Model stuck at a percentage | Check terminal for tracebacks. A common cause is Kp = −1 for an unsupported element |
| Depth mismatch warning | Verify the depth-unit dropdowns on both the Depth/Age and TE panels match your data |
| Pages/tabs not responding | Hard refresh; a JS error in the upload chain can freeze navigation |

---

## 7. Citation

If you use Dr Paleo, please cite:

> Hartland, A., Goswami, B., Park, J., Höpker, S.N., Torres Rojas, D., Liao, J., Fox, B.R.S., Marwan, N., Breitenbach, S.F.M. & Hu, C. (in review). Quantitative Holocene precipitation reconstruction from stalagmite trace metal kinetics reveals East Asian monsoon drivers. Preprint and DOI to follow.

Code and data: [github.com/waikatosci/paleodriprates](https://github.com/waikatosci/paleodriprates) · [doi:10.5281/zenodo.16392750](https://doi.org/10.5281/zenodo.16392750)

---

## 8. Contact

Adam Hartland — [adam.hartland@lincolnagritech.co.nz](mailto:adam.hartland@lincolnagritech.co.nz)

Lincoln Agritech Ltd / University of Waikato, New Zealand
