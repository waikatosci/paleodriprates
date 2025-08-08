PaleodripRates: Kinetic Proxy for Stalagmite Drip Rate and Precipitation Reconstruction

This repository contains code, data, and utilities for reconstructing cave drip rates and Holocene precipitation from stalagmite trace metals using a kinetic proxy model based on organic-metal complex dissociation. The methods are detailed in the manuscript "Quantitative Holocene Precipitation Reconstruction from Stalagmite Trace Metal Kinetics Reveals East Asian Monsoon Drivers" (Hartland et al., Nature Geoscience, 2025).
The proxy exploits the distributed dissociation kinetics of transition metals (e.g., Co, Ni) bound to organic ligands in dripwater, calibrated against modern data from Heshang Cave, China. Drip rates are inverted probabilistically via Bayesian methods, then chained to precipitation/temperature (P/T) regressions and temperature proxies for quantitative precipitation estimates (mm/yr) with uncertainties.
All code is written in Python and designed for reproducibility. Data and outputs are included where possible; full datasets are archived on Zenodo for persistent access.

Dependencies
Python 3.12+
Required libraries: numpy, pandas, scipy (integrate, optimize, interpolate, stats), matplotlib, openpyxl, Pillow (PIL), progressbar, bayprox (custom Bayesian library included in repo)
No additional installations are needed beyond standard pip-installable packages. For full execution, use the provided Excel data files.

Installation
Clone the repository:

git clone https://github.com/waikatosci/paleodriprates.git
cd paleodriprates

Install dependencies (if not already present):
pip install -r requirements.txt

(Note: bayprox is a custom module; it's included in the repo under lib/—add to PYTHONPATH if needed.)

Usage
1) Kinetic Model and Dissociation Kinetics:
Run model.py for numerical integration of dissociated fractions (e.g., log-normal distribution of rate constants).
Example: Compute labile fraction for given residence time (τ).

2) Bayesian Inversion for Drip Rates:
Use drip_rate_util.py to estimate drip rates from trace metal posteriors.
Input: Proxy records (Co/Ni concentrations) processed via BayProx for age-depth uncertainties.
Output: Drip rate PDFs (medians, percentiles) saved as drip_rate_percentiles.xlsx.

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
