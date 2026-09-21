# Extended-data figure scripts (Supplementary Figures 11–16)

Each folder carries a copy of the canonical native-resolution drip-rate
summary (`drip_rate_summary.csv`, 568 points, 16.66-target refit, 2026-09-21).
Depth-to-age conversion uses the HS4 chronology in
`dr_app/HS4_example_inputs/HS4_age_depth.csv` (30 ²³⁰Th ages plus the growth
surface at −51 yr BP) and, where a script asks for it, the interpolated 1-yr
age model written by the Dr Paleo run (`age_model.csv`, columns `depth`,
`age_yBP`).

Working invocations from each folder (paths relative to the repository root):

| Figure | Command |
|---|---|
| SF11 | `python Ex_D1.py --age_model <run>/age_model.csv` |
| SF12 | `python RQA_HS4_ensemble.py` (see `companion_analysis/RQA_HS4_ensemble.py`) |
| SF13 | `PYTHONPATH=../../manuscript_figures python ExD_3.py` (needs `age_model.csv`) |
| SF14 | `python ExD_4.py --tiepoints ../../dr_app/HS4_example_inputs/HS4_age_depth.csv` |
| SF15 | `python ExD_5.py --age_depth ../../dr_app/HS4_example_inputs/HS4_age_depth.csv --isotopes ../../Drip_rate.xlsx` (needs `astropy`) |
| SF16 | `python ExD_6.py --tiepoints ../../dr_app/HS4_example_inputs/HS4_age_depth.csv --age_model <run>/age_model.csv` |

`<run>` is the output folder of the canonical native run
(`manuscript_figures/external/run_logs/input_summary_native.csv` records its
settings). SF14 and SF15 reproduce the values quoted in their captions from
the files in this repository alone (SF15: drip-rate peaks at 371, 601, 774,
1,150, 1,641 and 3,816 yr; δ¹⁸O peaks at 353, 555, 746, 1,436 and 2,509 yr).
