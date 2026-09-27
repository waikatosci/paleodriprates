# v1.1.0 (28 September 2026): Nature Communications resubmission

State of the repository at resubmission of Hartland et al., "Decoupled infiltration and isotope
signals reveal a hidden East Asian monsoon megadrought" (NCOMMS-26-041445-T). Full detail of every
change since v1.0.2 is in PATCH_NOTES.md. In summary:

- Second-laboratory full ICP-MS suite added (`manuscript_figures/external/HS4_TE_full_suite_both_labs.xlsx`,
  with the measurement report released as Source Data); its offset from the primary run traced to the
  calcium matrix (`companion_analysis/crosslab_matrix_check.py`).
- The 5.2 ka event reported from Ni alone (`companion_analysis/run_single_metal.py`,
  `ni_co_agreement.py`): minimum 6.4 drips/min, ~70% below the mid-Holocene level; the joint Ni-Co
  minimum (2.8 drips/min, ~85%) kept as an upper bound.
- The second-laboratory run used as an independent replication of the event on Ni after correction
  for the Ca of each digest (`crosslab_replication.py`).
- Detrital, lithogenic and prior-carbonate-precipitation screens for the 5.2 and 8.2 ka intervals
  (`detrital_ternary_screen.py`, `events_differ_lithogenic.py`, `element_covariance.py`).
- Dripwater integration timescale and growth seasonality (`dripwater_integration_timescale.py`).
- One-step precipitation calibration on the Heshang P/T series (`manuscript_figures/build_precip_onestep.py`).
- Figures 5 and 6 as submitted (Fig. 5: dark median and interquartile lines, modern-baseflow star at the
  young end; Fig. 6 comparator named TraCE-21K-II).
- Documentation: README, `dr_app/QUICKSTART.md`, `dr_app/PRODUCTION_SETTINGS.md` (rev 4, 568-sample input),
  `manuscript_figures/external/README.md`, `supplementary_figures/README.md`, `CITATION.cff`, `.zenodo.json`.
- Withdrawn: the photographic composite of the polished section and everything that depended on it.

Not in this repository (deposited under 10.5281/zenodo.16392750): the age-propagated realisation
ensemble `drip_rate_realisations_ap.csv.gz`.
