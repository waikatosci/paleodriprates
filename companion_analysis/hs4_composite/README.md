# hs4_composite

Photographic composite of the polished HS4 section, used for Supplementary Figures 21–22 and Fig. 5c.

| script | does |
|---|---|
| `build_composite.py` | stitches the eight DSLR photographs (C. Hu; not in the repository, set `PHOTO_DIR`) onto the pencil depth scale: finds the scale line and its 5 cm ticks in each frame, fits depth to image row, resamples to 100 px/cm with the line vertical, aligns and colour-matches overlapping frames of the same piece, blends them, and cuts the 6 cm axial slice. Writes `manuscript_figures/external/HS4_composite.*` and `HS4_axial_slice.*` |
| `si_figures.py` | Supplementary Figure 21 (the composite, with the ²³⁰Th sample depths) and Supplementary Figure 22 (148–164 cm with primary-run and Ca-corrected second-run Ni and Co) |

Registration: tick residuals are below 1 mm in every frame (0.08–0.79 mm), and frames of the same piece agree to
0.1–0.5 mm where they overlap. The composite covers the stone from ~7 to 244 cm; the top 7 cm of the record and the base below 244 cm
were not photographed. Sample depths are plotted on the pencil scale as given.

The visible banding was also tested as a growth-rate measure (band spacing against the ²³⁰Th growth rate between
dates). At ~0.1 mm per pixel the photographs do not resolve annual layers (0.15–0.4 mm/yr), and the visible bands
(0.9–1.6 mm, 2–9 years each) do not track the ²³⁰Th growth rate (Spearman ρ −0.08 to 0.24, p ≥ 0.29), so no
lamina counts are reported.
