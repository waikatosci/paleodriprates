# `legacy/` — preserved historical files

This directory holds files that have been **superseded** in the active codebase but are kept here for **provenance**, reproducibility, and the historical record. Nothing in `legacy/` is imported or executed by the live code — everything still works without it. Treat this directory as a read-only archive.

If you are looking for the current canonical workflow, see the top-level `README.md` and `dr_app/QUICKSTART.md`.

---

## `legacy/scripts/` — superseded CLI drivers

| File | Replaced by | Why preserved |
|---|---|---|
| `Drip_rate.py` | `drip_rate_mc_realisations.py` (root) | Earliest CLI variant. Predates the parallel and full-realisations versions. |
| `Drip_rate_serial.py` | `drip_rate_mc_realisations.py` (root) | Single-core debug variant. Useful as a reference if anyone needs to step through the kinetic inversion serially without thread-pool noise. |
| `Drip_rate_parallel.py` | `drip_rate_mc_realisations.py` (root) | Parallel version that did **not** emit the full MC realisation ensemble (only percentile summaries). The `_fr` ("full realisations") successor — now `drip_rate_mc_realisations.py` — added the realisation export needed by `drip_rate_stationarity_tests.py`. |
| `RQA_HS4_ensemble.py` | `companion_analysis/RQA_HS4_ensemble.py` | Duplicate of the canonical file in `companion_analysis/` (one-line drift, almost certainly accidental). The companion-analysis copy is the live one. |

All four scripts read inputs from the root-level `Drip_rate.xlsx` and depend on the same shared core (`model.py`, `params.py`, `utils.py`, `drip_rate_util.py`). They will still run if invoked directly, but **they are not maintained** and may diverge from the live code over time.

## `legacy/env/` — obsolete environment files

| File | Status |
|---|---|
| `drip_rate_Windows7.yml` | Windows 7 reached end-of-life in January 2020. This conda env file targeted Win 7-era Python tooling and is preserved purely as a record of the historical build environment. The active env file is `drip_rate.yml` at the repo root. |

## `legacy/bayprox/` — superseded BayProX files

| File | Status |
|---|---|
| `proxyrecord_20170317.py` | Dated snapshot from 2017-03-17, predates the current `bayprox/proxyrecord.py`. Differs from the live version. Kept here in case the historical implementation is needed for reproducing pre-2018 results. |
| `Things.TODO` | Personal scratch / TODO file from early development. No code references it. Preserved because deleting personal notes from a public scientific repo can mask the messy reality of how the code was actually built. |

---

## When to update this directory

- **Add to it** if you're retiring a script but want to preserve provenance — the rule of thumb is "keep if it ever shipped, was cited, or might be needed to reproduce a published figure".
- **Don't delete from it** without strong reason. The cost of holding ~100 KB of old `.py` files is essentially zero; the cost of losing scientific provenance can be high.
- **Don't import from it** in active code. If you find yourself reaching into `legacy/`, that's a sign the file should be promoted back into active use (and updated/tested) rather than imported from the archive.

---

*Created 2026-04-09 as part of the `cleanup-2026-04` branch (T-0015).*
