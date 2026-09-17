# Supplementary Figures

`make_supp_figs_1_5.py` regenerates Supplementary Figures 1-5 from the released
data in house style (bare bold panel letters, no titles). Supplementary Figures
6-10 are produced by the companion and calibration scripts:

| Figure | Script |
|---|---|
| SF1-SF5 | `supplementary_figures/make_supp_figs_1_5.py` |
| SF6 | `companion_analysis/dripwater_source_stability.py` |
| SF7 | `companion_analysis/source_variation_propagation.py` |
| SF8 | `companion_analysis/detrital_ternary_screen.py` (two-panel enrichment summary, 2026-09-17) |
| SF9 | `companion_analysis/make_fig_events_differ.py` |
| SF10 | `calibration/make_fig_sigma_limits.py` (numerical demonstration of the sigma floor; supersedes `legacy/scripts/make_fig_censoring_mechanism.py`) |
| SF11-SF16 | `extended_data/Ex_Data_1..6` |

Outputs are written to `supplementary_figures/output/` (SF1-5) and
`manuscript_figures/output/` (SF6-10).
