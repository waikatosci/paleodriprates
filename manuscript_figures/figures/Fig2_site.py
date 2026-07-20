"""
Fig 2 — site map and drip-rate monitoring.

Reads from HS4_SourceData.xlsx:
  05_monitoring     drip-rate series
  Fig2_panelA_geo   site coordinates and monsoon-fringe control points

Panel A also uses two third-party gridded inputs, cited not embedded (see external/
and Fig2_panelB): GPCP v2.3 precipitation (precip.mon.mean.nc) and Natural Earth 10m
relief (SR_LR.tif). The 12 plotted monthly rainfall medians are archived on
Fig2_panelB and cross-checked against the netCDF at run time.
"""

import os as _os, sys as _sys
# Resolve everything against the repo root so generate_figures.py can run each script
# from one place, and so a figure run from its own directory behaves identically.
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_sys.path.insert(0, _ROOT)
_os.chdir(_ROOT)
_os.makedirs('output', exist_ok=True)

import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LongitudeFormatter, LatitudeFormatter
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import matplotlib.image as mpimg

matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans']
matplotlib.rcParams['font.size'] = 7

# ── THE SINGLE SOURCE ────────────────────────────────────────────────────────
SOURCE = 'HS4_SourceData.xlsx'          # canonical workbook
DATA   = 'external'                        # cited gridded products only (see header)
OUT_STEM = 'output/Fig2_site'

# ── AUTHOR-DIRECTED DISPLAY SETTINGS (2026-07-15; unchanged) ─────────────────
LEVELS_MAX   = 12.0
RELIEF_ALPHA = 0.55
RELIEF_VMIN  = 120
RELIEF_VMAX  = 206

COL_MODERN  = '#0D47A1'
COL_HOL     = '#E65100'
COL_LGM     = '#00ACC1'
COL_SITES   = 'seagreen'
COL_HESHANG = 'maroon'

# ── read the workbook ────────────────────────────────────────────────────────
def _find_header(sheet, first_col_name):
    """Locate a sheet's header row by name rather than by a hardcoded skiprows.

    Sheets grow notes and provenance blocks above the table as findings accumulate —
    that is the workbook doing its job. Binding to a fixed row number makes every
    script hostage to that; binding to the column name does not. (Learned the hard
    way: a note added to 05_monitoring for Fig 4 broke this script's skiprows=4.)
    """
    raw = pd.read_excel(SOURCE, sheet_name=sheet, header=None)
    hits = raw.index[raw[0].astype(str).str.strip() == first_col_name]
    if len(hits) == 0:
        raise ValueError(f"header '{first_col_name}' not found in sheet '{sheet}'")
    return int(hits[0])


# 05_monitoring is now built on the MASTER (HS4_Drip_discharge.xlsx, mL/min), which
# arrived with the Fig 6 bundle. Fig 2b plots the COUNTED drips/min series, which is a
# different quantity from the master's volumetric mL/min and is carried alongside it in
# 'counted_drips_per_min'. Reading that column keeps this figure byte-identical to the
# shipped one. (Sourcing Fig 2b from the master instead would move the July peak
# 22.8 -> 30.0 drips/min — a real change to a published panel, and an author decision.)
mon = pd.read_excel(SOURCE, sheet_name='05_monitoring',
                    skiprows=_find_header('05_monitoring', 'date'))
mon = mon.dropna(subset=['counted_drips_per_min'])
df = pd.DataFrame({'date': pd.to_datetime(mon['date']),
                   'DR': mon['counted_drips_per_min'].astype(float)})
df['month'] = df['date'].dt.month
median_dr = df.groupby('month')['DR'].median()
print(f'[source] 05_monitoring: {len(df)} drip points, '
      f'{df.date.min():%Y-%m-%d} to {df.date.max():%Y-%m-%d}')

_pb = pd.read_excel(SOURCE, sheet_name='Fig2_panelB',
                    skiprows=_find_header('Fig2_panelB', 'month')).dropna(subset=['month'])
_pb = _pb[_pb['month'].apply(lambda v: isinstance(v, (int, float)))].head(12)
median_precip_book = pd.Series(_pb['median_rainfall_mm_day'].astype(float).values,
                               index=_pb['month'].astype(int).values)

geo = pd.read_excel(SOURCE, sheet_name='Fig2_panelA_geo',
                    skiprows=_find_header('Fig2_panelA_geo', 'site'))
geo = geo.dropna(subset=['site'])
sites = {r['site']: (float(r['lon_E']), float(r['lat_N']))
         for _, r in geo[geo['type'] == 'cave'].iterrows()}
marine_cores = {r['site']: (float(r['lon_E']), float(r['lat_N']))
                for _, r in geo[geo['type'] == 'marine core'].iterrows()}
_lk = geo[geo['type'] == 'lake'].iloc[0]
lake_yoa = (float(_lk['lon_E']), float(_lk['lat_N']))
print(f'[source] Fig2_panelA_geo: {len(sites)} caves, {len(marine_cores)} cores, 1 lake')

fr = pd.read_excel(SOURCE, sheet_name='Fig2_panelA_geo',
                   skiprows=_find_header('Fig2_panelA_geo', 'lon_E')).head(7)
fringe_lon_points = fr['lon_E'].astype(float).values
modern_lat_mean   = fr['modern_lat_N'].astype(float).values
holocene_shift    = fr['holocene_shift_deg'].astype(float).values
lgm_shift         = fr['lgm_shift_deg'].astype(float).values
print(f'[source] fringe control points: {len(fringe_lon_points)}')

# ── GPCP precipitation (cited gridded product; unchanged logic) ──────────────
ds = xr.open_dataset(f'{DATA}/precip.mon.mean.nc', engine='h5netcdf')
precip = ds['precip'].sel(time=slice('1998-01-01', '2025-07-01')).mean(dim='time')
# prime-meridian seam fix (0-360 -> -180..180), as parent
precip = precip.assign_coords(lon=(((precip.lon + 180) % 360) - 180)).sortby('lon')
fine_lat = np.linspace(precip.lat.min().item(), precip.lat.max().item(), num=precip.lat.size * 2)
precip_interp = precip.interp(lat=fine_lat, method='linear')

heshang_lat, heshang_lon = 30.45, 110.417
precip_ts = ds['precip'].sel(time=slice('1998-01-01', '2025-07-01')).sel(
    lat=heshang_lat, lon=heshang_lon, method='nearest')
median_precip = precip_ts.groupby('time.month').median()

# the archived copy must agree with the live derivation
_dev = np.abs(median_precip.values - median_precip_book.reindex(median_precip.month.values).values)
assert _dev.max() < 1e-5, f'workbook rainfall disagrees with netCDF (max dev {_dev.max():.2e})'
print(f'[check] Fig2_panelB rainfall matches netCDF (max dev {_dev.max():.2e})')

# ── figure / map ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(7, 3))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-20, 140, 0, 45], crs=ccrs.PlateCarree())
fig.subplots_adjust(left=0.05, right=0.7, bottom=0.25, top=0.95)

img = mpimg.imread(f'{DATA}/SR_LR.tif')
ny, nx = img.shape[:2]
lon0, lon1, lat0, lat1 = -22, 142, -2, 47
c0, c1 = int((lon0 + 180) / 360 * nx), int((lon1 + 180) / 360 * nx)
r0, r1 = int((90 - lat1) / 180 * ny), int((90 - lat0) / 180 * ny)
ax.imshow(img[r0:r1:2, c0:c1:2], extent=[lon0, lon1, lat0, lat1], origin='upper',
          transform=ccrs.PlateCarree(), alpha=RELIEF_ALPHA, cmap='gray',
          vmin=RELIEF_VMIN, vmax=RELIEF_VMAX)

cmap = plt.cm.Blues
levels = np.arange(0, LEVELS_MAX + 0.25, 0.25)
cf = precip_interp.plot.contourf(ax=ax, transform=ccrs.PlateCarree(), levels=levels,
                                 cmap=cmap, extend='max', alpha=0.5, add_colorbar=False)
cbar = plt.colorbar(cf, ax=ax, orientation='horizontal', pad=0.1, aspect=30, shrink=0.8)
cbar.set_label('Mean Precipitation (mm/day, GPCP 1998-2025)')

ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
gl = ax.gridlines(draw_labels=True, linestyle='--', alpha=0.3, linewidth=0.5)
gl.top_labels = False
gl.right_labels = False
gl.xformatter = LongitudeFormatter()
gl.yformatter = LatitudeFormatter()

# ── monsoon fringes ──────────────────────────────────────────────────────────
fringe_lon = np.linspace(-20, 150, 200)
modern_lat = interp1d(fringe_lon_points, modern_lat_mean, kind='cubic')(fringe_lon)
holocene_lat = interp1d(fringe_lon_points, modern_lat_mean + holocene_shift, kind='cubic')(fringe_lon)
lgm_lat = interp1d(fringe_lon_points, modern_lat_mean + lgm_shift, kind='cubic')(fringe_lon)
band_width_modern, band_width_paleo = 2.5, 3.5


def plot_gradient_band(ax, x, y_mean, band_width, color, alpha_max, steps=10):
    for i in range(steps):
        frac = (steps - i) / steps
        ax.fill_between(x, y_mean - band_width * frac, y_mean + band_width * frac,
                        color=color, alpha=alpha_max / steps, linewidth=0)


plot_gradient_band(ax, fringe_lon, lgm_lat, band_width_paleo, COL_LGM, alpha_max=0.1)
plot_gradient_band(ax, fringe_lon, holocene_lat, band_width_paleo, COL_HOL, alpha_max=0.1)
plot_gradient_band(ax, fringe_lon, modern_lat, band_width_modern, COL_MODERN, alpha_max=0.3)

lgm_line, = ax.plot(fringe_lon, lgm_lat, color=COL_LGM, linewidth=0.5, linestyle=':',
                    label='LGM Monsoon Fringe (±3.5°)')
holocene_line, = ax.plot(fringe_lon, holocene_lat, color=COL_HOL, linewidth=0.5, linestyle='--',
                         label='Holocene Monsoon Fringe (±3.5°)')
modern_line, = ax.plot(fringe_lon, modern_lat, color=COL_MODERN, linewidth=0.75, linestyle='-',
                       label='Modern Monsoon Fringe (±2.5°)')

# ── sites ────────────────────────────────────────────────────────────────────
heshang_handle = None
for site, coord in sites.items():
    size = 8 if site == 'Heshang Cave' else 5
    if site == 'Heshang Cave':
        heshang_handle = ax.scatter(coord[0], coord[1], marker='*', s=size ** 2,
                                    color=COL_HESHANG, label=f'{site}', zorder=11)
    else:
        ax.scatter(coord[0], coord[1], marker='*', s=size ** 2, facecolor='none',
                   edgecolor=COL_SITES, linewidths=0.5, zorder=10)
for core, coord in marine_cores.items():
    ax.scatter(coord[0], coord[1], marker='s', s=25, facecolor='none',
               edgecolor=COL_SITES, linewidths=0.5, zorder=10)
ax.scatter(lake_yoa[0], lake_yoa[1], marker='o', s=25, facecolor='none',
           edgecolor=COL_SITES, linewidths=0.5, zorder=10)

initials = {'Soreq Cave': 'SC', 'Tham Doun Mai Cave': 'TDMC', 'Dongge Cave': 'DC',
            'Sanbao Cave': 'SBC', 'Hulu Cave': 'HC', 'Qunf Cave': 'QC',
            'Jeita Cave': 'JC', 'Jiuxian Cave': 'JXC', 'Mawmluh Cave': 'MC'}
marine_initials = {'Hole 658C': 'H658C', 'GeoB9508-5': 'GB9508'}
dx, dy = 1.5, 0.3
label_font = dict(fontsize=5, color='black', fontweight='light', ha='left', va='center', zorder=12)
for site, coord in sites.items():
    if site != 'Heshang Cave':
        ax.text(coord[0] + dx, coord[1] + dy, initials[site], **label_font)
for core, coord in marine_cores.items():
    ax.text(coord[0] + dx, coord[1] + dy, marine_initials[core], **label_font)
ax.text(lake_yoa[0] + dx, lake_yoa[1] + dy, 'LY', **label_font)

ax.legend(handles=[lgm_line, holocene_line, modern_line, heshang_handle],
          bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=6, borderaxespad=0.)

banner_text = (
    "Cave sites (open stars): SC = Soreq Cave, TDMC = Tham Doun Mai Cave, DC = Dongge Cave, "
    "SBC = Sanbao Cave, HC = Hulu Cave, QC = Qunf Cave, JC = Jeita Cave\n"
    "Marine cores (open squares): H658C = Hole 658C, GB9508 = GeoB9508-5\n"
    "Lake sites (open circle): LY = Lake Yoa")
fig.text(0.05, 0.05, banner_text, fontsize=6, va='bottom', ha='left', wrap=True)
ax.text(-0.05, 1.05, 'A', transform=ax.transAxes, fontsize=10, fontweight='bold', va='top')

# ── Panel B ──────────────────────────────────────────────────────────────────
ax_inset = fig.add_axes([0.76, 0.22, 0.14, 0.42])
months = np.arange(1, 13)
months_fine = np.linspace(1, 12, 200)

ax2 = ax_inset.twinx()
f_precip = interp1d(months, median_precip.values, kind='linear')
ax2.fill_between(months_fine, 0, f_precip(months_fine), color='steelblue', alpha=0.2, zorder=1)
ax2.plot(months_fine, f_precip(months_fine), color='steelblue', label='Median rainfall', zorder=2)

f_dr = interp1d(months, median_dr.values, kind='linear')
ax_inset.plot(months_fine, f_dr(months_fine), color='black', linestyle='-',
              label='Median drip rate', zorder=3)
ax_inset.plot(months, median_dr.values, color='black', marker='o', linestyle='none',
              markersize=3, zorder=4)

assert f_dr(months_fine).max() <= median_dr.values.max() + 1e-9, 'drip interp exceeds observed!'
assert f_precip(months_fine).max() <= median_precip.values.max() + 1e-9, 'precip interp exceeds observed!'
print(f'Panel B check: drip interp max {f_dr(months_fine).max():.2f} <= observed max '
      f'{median_dr.values.max():.2f}; rainfall interp max {f_precip(months_fine).max():.2f} '
      f'<= observed max {median_precip.values.max():.2f}')

ax_inset.set_ylim(0, 26)
ax2.set_ylim(0, 8)
ax_inset.set_xlim(0.5, 12.5)
ax_inset.set_xticks(months)
ax_inset.set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'], fontsize=6)
ax_inset.set_ylabel('Drip rate (drips/min)', color='black', fontsize=6)
ax2.set_ylabel('Rainfall (mm/day)', color='steelblue', fontsize=6)
ax_inset.tick_params(axis='y', colors='black', labelsize=6)
ax2.tick_params(axis='y', colors='steelblue', labelsize=6)
ax_inset.tick_params(axis='x', labelsize=6)
ax_inset.text(-0.35, 1.05, 'B', transform=ax_inset.transAxes, fontsize=10,
              fontweight='bold', va='top')

for fmt in ('pdf', 'png'):
    fig.savefig(f'{OUT_STEM}.{fmt}', dpi=300, bbox_inches='tight')
    print(f'Saved -> {OUT_STEM}.{fmt}')
plt.close(fig)
print('Done.')
