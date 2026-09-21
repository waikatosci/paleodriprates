# Round-2 revision patch — second-laboratory full ICP-MS suite

Add-on to the paleodriprates repository at commit `fb3fc80`.

## Files to add (new)

| path | size | sha256 |
|---|---|---|
| `manuscript_figures/external/HS4_TE_full_suite_both_labs.xlsx` | 51,375 B | `72dd17205b3fd0b133186e41dad8cb213a4c7d18003ff2fb57f79b061750b1e7` |
| `manuscript_figures/external/HS4_TE_full_suite_both_labs.README.md` | 2,476 B | (this note) |
| `companion_analysis/crosslab_matrix_check.py` | 12,508 B | (this note) |

## Files to update (edit in place)

- `calibration/excluded_points.csv` — reason field on the 20 second-lab rows
  updated from "certified-standard failure" to the Ca-matrix explanation
  (the corrected file is `excluded_points.csv` in the delivery bundle).
- `manuscript_figures/external/README.md` — add one row under the
  "Added for the revision" table:

  ```
  | `HS4_TE_full_suite_both_labs.xlsx` | `companion_analysis/crosslab_matrix_check.py`, Supp. Methods 14.3–14.4, Supp. Fig. 18 | full trace-element suite (Li–U) for 21 HS4 powders and four silicate reference materials; the 2019 second-laboratory run that supplied the twenty samples now withdrawn from the input, plus the paired 2009 primary-run values at the three shared depths (see companion README) |
  ```

- `companion_analysis/crosslab_sensitivity.py` — replace the docstring
  reference to a "certified-standard failure" with a pointer to
  `crosslab_matrix_check.py` and the Ca-matrix interference. The
  sensitivity analysis and its numbers do not change.

## Commit message

```
Add the second-laboratory full ICP-MS suite and identify the Ca matrix as the cause of its offset

The second laboratory (Wuhan, November 2019) measured the full trace-element
suite for 21 HS4 powders (12 across the 5.2 ka event, 155.2-157.8 cm; 8 near
the growth surface, 7.0-8.5 cm) and returned AGV-2, BHVO-2 and BCR-2 within a
few percent of certified values. It nonetheless offsets from the primary run
at the three depths analysed by both laboratories. The cause is the calcium
matrix: outside the 5.2 ka core, Co, Ni and Fe in the second run are linear
functions of Ca in the digest (r = 0.98, 0.96, 0.97) with a contribution at
calcite Ca of ~1.6 ppm Co, ~5.8 ppm Ni and ~6,900 ppm Fe (~0.7%). The Fe is
the 40Ca16O+ polyatomic interference at m/z 56; the Co and Ni baselines are
consistent with the corresponding Ca-bearing ions at m/z 59 and 60. The
primary run's standards were matrix-matched to high-purity Ca and do not
carry it. The 20 samples remain withdrawn from the inversion for that
reason.

The same file also provides Al, Th, Ti, Zr and Mn at the 5.2 ka horizon,
which the primary run did not measure. Al up to 38 ppm, Th up to 0.013 ppm
and Ti up to 7 ppm bound any detrital contribution at below 1 percent of the
observed Co and Ni excesses; the excess Co and Ni scale with Mn (Co/Mn 0.022,
Ni/Mn 0.039) with no lithogenic response, the signature of a residence-time
excursion in four organically bound metals. This is the mass balance quoted
in Supplementary Methods 14.4 and Supplementary Figure 18.

  add:   manuscript_figures/external/HS4_TE_full_suite_both_labs.xlsx
  add:   manuscript_figures/external/HS4_TE_full_suite_both_labs.README.md
  add:   companion_analysis/crosslab_matrix_check.py
  edit:  manuscript_figures/external/README.md
  edit:  calibration/excluded_points.csv
  edit:  companion_analysis/crosslab_sensitivity.py    # docstring only

Reproduces every value in Supplementary Methods 14.3, 14.4 and Supplementary
Figure 18 by running `python companion_analysis/crosslab_matrix_check.py`.

Data provided by C. Hu, 21 September 2026.
```
