"""
Fig 7 — HS4 Co/Ni drip-rate × IRM_soft correspondence, with IRM event detection.

Reads:
  external/pdf_heatmap_hr.json            drip-rate PDF heatmap (σ = π/√6)
  external/drip_rate_summary_{hr,lr}.csv  drip-rate percentiles
  external/HS4_Zhu2017_IRMsoft_flux.csv   IRM_soft flux record (Zhu 2017)
  HS4_SourceData.xlsx[02_chronology]      age model and U-Th tie-points

Ages: Fig 7 uses np.interp (clamps), placing the 255.50 cm sample at 9.446 ka
(Fig 5 extrapolates to 9.604 ka — an author decision recorded on 02_chronology).

Z-score event detection on IRM_soft flux (Mayewski et al. 2004 approach):
      * Interpolate IRM onto a regular 50-yr grid
      * Compute a centred ±750-yr rolling mean + SD
      * Flag z > 2 exceedances; merge contiguous flags (gap tolerance 200 yr)
      * Draw vertical brown bars (behind both series) at detected events
  - Detection is IRM-only. IRM_soft tracks flood-energy / extreme-event delivery,
    modulated by antecedent conditions (which can be dry OR wet).
  - Results at ±750-yr / z>2 — three events:
      5.24–5.29 ka (z=2.80): drip-rate minimum; drought-primed flood energy.
      4.14 ka       (z=2.25): coincides with the middle-Yangtze WET expression of
                              the 4.2 ka event (Meghalayan) — regional pluvial /
                              paleoflood interval, NOT drought at Heshang.
      0.99 ka       (z=2.20): late-Holocene drip minimum; drought-primed.
    The common control is high-magnitude episodic transport, not drying per se.

References:
  Mayewski et al. (2004) Quat Sci Rev 23:1257–1278 — rolling z-score for
  Holocene rapid climate events.
  Istanbulluoglu & Bras (2006) WRR 42 W06418 — ecohydrological mechanism.
  Hartland et al. (2012) Chem Geol 304–305:68 — high/low flux transport modes.
"""

import os as _os, sys as _sys
# Resolve everything against the repo root so generate_figures.py can run each script
# from one place, and so a figure run from its own directory behaves identically.
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_sys.path.insert(0, _ROOT)
_os.chdir(_ROOT)
_os.makedirs('output', exist_ok=True)

import json, numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FuncFormatter
from matplotlib.transforms import blended_transform_factory
import ngeo_style as ng
ng.apply_style()

# ── the single source ──────────────────────────────────────────────────
SOURCE = "HS4_SourceData.xlsx"
HM     = "external/pdf_heatmap_hr.json"
IRM_FILE = "external/HS4_Zhu2017_IRMsoft_flux.csv"


def _find_header(sheet, first_col_name):
    """Locate a sheet's header row by name rather than a hardcoded skiprows."""
    raw = pd.read_excel(SOURCE, sheet_name=sheet, header=None)
    hits = raw.index[raw[0].astype(str).str.strip() == first_col_name]
    if len(hits) == 0:
        raise ValueError(f"header '{first_col_name}' not found in '{sheet}'")
    return int(hits[0])


def _sheet(name, first_col):
    return pd.read_excel(SOURCE, sheet_name=name,
                         skiprows=_find_header(name, first_col)).dropna(how="all")


# ── age model ← 02_chronology (was age_model.csv) ──────────────────────
am = _sheet("02_chronology", "depth_cm").rename(columns={"depth_cm": "depth"})
am = am[pd.to_numeric(am["depth"], errors="coerce").notna()]

# np.interp CLAMPS beyond the deepest date. That is why Fig 7 places the 255.50 cm
# sample at 9.4460 ka while Fig 5, which extrapolates, places it at 9.6035 ka.
# The 255.50 cm age is an author decision recorded on 02_chronology.
def d2a(depth): return np.interp(depth, am['depth'], am['age_yBP']) / 1000.

# ── heatmap data ────────────────────────────────────────────────────────
with open(HM) as f: hm = json.load(f)
V_span    = np.array(hm['V_span'])
hm_depths = np.array(hm['ages'])
V_pdf     = np.array(hm['V_pdf']).T        # (n_age, n_drip) after transpose
hm_ages   = d2a(hm_depths)

ok = np.isfinite(hm_ages) & (hm_depths > 0)
hm_ages, V_pdf = hm_ages[ok], V_pdf[ok]
for i in range(V_pdf.shape[0]):            # normalise per age slice
    mx = V_pdf[i].max()
    if mx > 0: V_pdf[i] /= mx
idx = np.argsort(hm_ages)
hm_ages, V_pdf = hm_ages[idx], V_pdf[idx]

ae = np.zeros(len(hm_ages) + 1)
ae[1:-1] = 0.5 * (hm_ages[:-1] + hm_ages[1:])
ae[0]    = hm_ages[0]  - (hm_ages[1]  - hm_ages[0])  / 2
ae[-1]   = hm_ages[-1] + (hm_ages[-1] - hm_ages[-2]) / 2
de = np.zeros(len(V_span) + 1)
de[1:-1] = 0.5 * (V_span[:-1] + V_span[1:])
de[0]    = max(0, V_span[0] - (V_span[1] - V_span[0]) / 2)
de[-1]   = V_span[-1] + (V_span[-1] - V_span[-2]) / 2
AGE_E, DRIP_E = np.meshgrid(ae, de)

# ── drip rate ──────────────────────────────────────────────────────────
dhr = pd.read_csv("external/drip_rate_summary_hr.csv")
dhr.columns = [c.lower() for c in dhr.columns]
dhr["ka"] = d2a(dhr["depth"].values)
dhr = dhr[np.isfinite(dhr["ka"]) & (dhr["pc50"] > 0)].sort_values("ka")

dlr = pd.read_csv("external/drip_rate_summary_lr.csv")
dlr.columns = [c.lower() for c in dlr.columns]
dlr["ka"] = d2a(dlr["depth"].values)
dlr = dlr[np.isfinite(dlr["ka"]) & (dlr["pc50"] > 0)].sort_values("ka")

# ── IRM ────────────────────────────────────────────────────────────────
irm = pd.read_csv(IRM_FILE).sort_values("age_mid_kaBP")
irm["flux9"] = irm["IRMsoft_flux_Am2_per_yr"] * 1e9

# ── IRM z-score event detection (Mayewski et al. 2004 approach) ────────
# Interpolate onto regular 50-yr grid; centred ±500-yr rolling window.
GRID_STEP  = 0.05   # ka
HALF_WIN   = 0.75   # ka  → 1500-yr total window
Z_THRESH   = 2.0
GAP_TOL    = 0.20   # ka — merge flags closer than this

grid  = np.arange(irm.age_mid_kaBP.min(),
                  irm.age_mid_kaBP.max() + GRID_STEP, GRID_STEP)
f9g   = np.interp(grid, irm.age_mid_kaBP.values, irm.flux9.values)
s     = pd.Series(f9g, index=grid)
half  = int(HALF_WIN / GRID_STEP)
rol_mu = s.rolling(2*half+1, center=True, min_periods=5).mean()
rol_sd = s.rolling(2*half+1, center=True, min_periods=5).std()
z_ser  = (s - rol_mu) / rol_sd

exceed = z_ser[z_ser > Z_THRESH].index.values
event_windows = []
if len(exceed):
    gaps   = np.diff(exceed)
    splits = np.where(gaps > GAP_TOL)[0] + 1
    segs   = np.split(exceed, splits)
    for seg in segs:
        if len(seg):
            lo, hi = seg[0], seg[-1]
            pz = z_ser[(z_ser.index >= lo) & (z_ser.index <= hi)].max()
            event_windows.append((lo, hi, pz))

print(f"IRM z-score events (win=±{HALF_WIN:.2f}ka, z>{Z_THRESH}):")
for lo, hi, pz in event_windows:
    print(f"  {lo:.2f}–{hi:.2f} ka  peak z={pz:.2f}")

# ── U-Th tie-points ─────────────────────────────────────────────────────
# tie-points ← 02_chronology (was HS4_UTh_tiepoints.csv), so that — as in Fig 5 —
# the plotted tie-points and the age model cannot diverge.
uth = _sheet("02_chronology", "depth_cm").rename(
    columns={"age_yBP": "age_bp", "err2s_yr": "err2s"})
uth = uth[pd.to_numeric(uth["age_bp"], errors="coerce").notna()]
uth["ka"]     = uth["age_bp"]  / 1000.
uth["err_ka"] = uth["err2s"]   / 1000.
uth = uth[(uth.ka >= 0) & (uth.ka <= 9.5) & (uth.err2s > 1)]

# ── helpers ─────────────────────────────────────────────────────────────
def gsm(x, y, g, s):
    return np.array([
        np.sum(np.exp(-0.5*((x-t)/s)**2)*y) /
        np.sum(np.exp(-0.5*((x-t)/s)**2)) for t in g])

def irm_to_drip(v):
    """Map IRM flux (×10⁻⁹) to primary drip-rate axis coordinates.
    IRM ylim (bot=6, top=-2) ↔ drip ylim (bot=0, top=40).
    """
    return (np.asarray(v) - (-2.0)) / (6.0 - (-2.0)) * (0.0 - 40.0) + 40.0

XLIM     = (9.5, 0.0)
BASEFLOW = 17.0
COL_MEDIAN = ng.COL_BG_900
IRM_COL    = "#FF8F00"
IRM_RING   = "#B35900"
EVT_COL    = "#795548"   # brown — behind both series

# ── figure ──────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(180/25.4, 85/25.4))
ax  = fig.add_subplot(111)
fig.subplots_adjust(left=0.16, right=0.82, top=0.90, bottom=0.22)

# ── IRM event bars (lowest zorder — behind everything) ─────────────────
for lo, hi, pz in event_windows:
    # minimum half-width: 100 yr so narrow events remain visible
    hw = max((hi - lo) / 2, 0.15)
    ctr = (lo + hi) / 2
    ax.axvspan(ctr - hw, ctr + hw,
               color=EVT_COL, alpha=0.45, lw=0, zorder=0.3)

# ── heatmap ─────────────────────────────────────────────────────────────
hm_cmap = LinearSegmentedColormap.from_list('ngeo_fig2blue_light', [
    (1.0,  1.0,  1.0,  0.00),
    (0.97, 0.98, 1.00, 0.15),
    (0.84, 0.90, 0.96, 0.28),
    (0.51, 0.73, 0.86, 0.45),
    (0.20, 0.51, 0.75, 0.60),
    (0.05, 0.34, 0.63, 0.72),
    (0.03, 0.19, 0.42, 0.82),
], N=256)
ax.pcolormesh(AGE_E, DRIP_E, V_pdf.T,
              cmap=hm_cmap, vmin=0, vmax=1,
              shading='flat', rasterized=True, zorder=1)

# ── drip rate (medians only, no IQR lines) ──────────────────────────────
ax.plot(dhr.ka, dhr.pc50, color=COL_MEDIAN,
        lw=0.3, alpha=0.25, zorder=3)
l_drip, = ax.plot(dlr.ka, dlr.pc50, color=COL_MEDIAN,
                  lw=0.5, zorder=4,
                  label="Co/Ni drip rate (this study)")
ax.axhline(BASEFLOW, color=ng.COL_BG_400, lw=0.5,
           ls=(0,(4,3)), zorder=1.5)
ax.set_ylim(0, 40)
ax.set_ylabel("Drip rate (drips min$^{-1}$)")
ax.tick_params(axis='y')
ax.set_xlabel("Age (ka BP)")

# ── IRM secondary axis (labels only) ───────────────────────────────────
ax2 = ax.twinx()
ax2.set_ylim(6.0, -2.0)
ax2.yaxis.set_major_formatter(FuncFormatter(
    lambda v, _: f"{v:.0f}" if v >= 0 else ""))
ax2.set_yticks([0, 1, 2, 3, 4, 5, 6])
ax2.set_ylabel("IRM$_{soft}$ flux (10$^{-9}$ A m$^2$ yr$^{-1}$)",
               color=IRM_COL)
ax2.tick_params(axis='y', colors=IRM_COL)

# ── IRM series on primary axis at top zorder ────────────────────────────
gi    = np.arange(0.15, 8.5, 0.05)
IRM_Z = 10
ax.plot(irm.age_mid_kaBP, irm_to_drip(irm.flux9),
        color=IRM_COL, lw=0.5, alpha=0.55, zorder=IRM_Z)
ax.scatter(irm.age_mid_kaBP, irm_to_drip(irm.flux9),
           s=9, color=IRM_COL, edgecolors=IRM_RING,
           linewidths=0.6, alpha=0.9, zorder=IRM_Z + 0.5)
l_irm, = ax.plot(gi,
                 irm_to_drip(gsm(irm.age_mid_kaBP.values,
                                 irm.flux9.values, gi, 0.6)),
                 color=IRM_COL, lw=1.8, zorder=IRM_Z + 1,
                 label="IRM$_{soft}$ flux (Zhu 2017)")

# ── shared styling ──────────────────────────────────────────────────────
ng.draw_events([ax], ng.EVENTS_FIG6)
for a in (ax, ax2): a.set_xlim(*XLIM)
ax.axvline(5.19, color=ng.COL_BG_600, lw=0.4,
           ls=(0,(3,2)), alpha=0.6, zorder=0.7)
ng.style_ax(ax); ng.style_ax(ax2)
ax.patch.set_visible(False)

# ── U-Th tie-points + errorbars above the top spine ────────────────────
trans    = blended_transform_factory(ax.transData, ax.transAxes)
MARKER_Y = 1.04
LABEL_Y  = 1.09
ax.scatter(uth.ka, np.full(len(uth), MARKER_Y), marker="v", s=5,
           color=ng.COL_BG_700, zorder=8, clip_on=False, transform=trans)
ax.errorbar(uth.ka, np.full(len(uth), MARKER_Y), xerr=uth.err_ka,
            fmt="none", ecolor=ng.COL_BG_500, elinewidth=0.4,
            capsize=1.2, capthick=0.4, zorder=7, clip_on=False,
            transform=trans)
ax.text(9.35, LABEL_Y, "²³⁰Th (±2σ)", fontsize=4.8,
        color=ng.COL_BG_600, va="center", ha="left",
        clip_on=False, transform=trans)

# ── legend ──────────────────────────────────────────────────────────────
leg = ax.legend(handles=[l_drip, l_irm],
                loc="lower left", frameon=False, fontsize=5.0,
                handlelength=1.4, handletextpad=0.4,
                labelspacing=0.3, borderaxespad=0.6)

fig.canvas.draw()
inv  = ax.transData.inverted()
bbox = leg.get_window_extent(renderer=fig.canvas.get_renderer())
x0d, _ = inv.transform((bbox.x0, bbox.y0))
x1d, _ = inv.transform((bbox.x1, bbox.y1))
right_ka = min(x0d, x1d)
print(f"Legend right edge: {right_ka:.2f} ka  "
      f"({'OK' if right_ka >= 6 else 'WARNING: crosses 6 ka'})")

fig.savefig("output/Fig7_events.pdf", dpi=600, bbox_inches="tight",
            metadata={"CreationDate": None})
fig.savefig("output/Fig7_events.png", dpi=600, bbox_inches="tight")
print("saved output/Fig7_events.pdf/png")
