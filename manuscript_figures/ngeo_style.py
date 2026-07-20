"""
ngeo_style.py
=============
Shared Nature Communications house style for Hartland et al.
Import this module in all figure scripts:

    from ngeo_style import *

Provides: rcParams setup, colour palette, style_ax(), draw_events().
"""
import matplotlib.pyplot as plt
import numpy as np

# ══════════════════════════════════════════════════════════════════════
# rcParams — Nature Communications house style
# ══════════════════════════════════════════════════════════════════════
NGEO_RC = {
    'font.family':        'sans-serif',
    'font.sans-serif':    ['Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size':          6.5,
    'axes.labelsize':     6.5,
    'axes.labelweight':   'normal',
    'axes.titlesize':     6.5,
    'axes.linewidth':     0.25,
    'xtick.labelsize':    5.5,
    'ytick.labelsize':    5.5,
    'xtick.major.width':  0.25,
    'ytick.major.width':  0.25,
    'xtick.major.size':   2.5,
    'ytick.major.size':   2.5,
    'xtick.minor.width':  0.2,
    'ytick.minor.width':  0.2,
    'xtick.minor.size':   1.2,
    'ytick.minor.size':   1.2,
    'xtick.direction':    'in',
    'ytick.direction':    'in',
    'legend.fontsize':    5.5,
    'figure.dpi':         300,
    'pdf.fonttype':       42,
    'ps.fonttype':        42,
    'savefig.dpi':        300,
    'savefig.bbox':       'tight',
}

def apply_style():
    """Apply NGeo rcParams. Call once at script start."""
    plt.style.use('default')
    plt.rcParams.update(NGEO_RC)

# ══════════════════════════════════════════════════════════════════════
# COLOUR PALETTE — 3-hue system
#   Blue/indigo family: HS4 data / Ni / kinetic proxy
#   Brown family: comparisons / Co secondary / δ¹⁸O differencing
#   Charcoal:     reference / models / axes
# ══════════════════════════════════════════════════════════════════════

# ── Primary data (HS4 / Ni / kinetic proxy — matched to Fig 2 blue) ──
COL_TEAL_900  = '#072F6B'   # Fig 2 darkest navy
COL_TEAL_800  = '#0C57A0'   # Fig 2 dark blue — primary line colour
COL_TEAL_600  = '#3282BE'   # Fig 2 mid blue
COL_TEAL_400  = '#4F9BCB'   # Fig 2 steel blue
COL_TEAL_300  = '#82BADB'   # Fig 2 sky blue
COL_TEAL_200  = '#D6E6F4'   # Fig 2 pale blue — IQR / CI fill
COL_TEAL_50   = '#F7FAFF'   # Fig 2 ice white — 90% CI fill

# ── Secondary data (Co / δ¹⁸O / brown family) ──
COL_BROWN_500 = '#795548'   # δ¹⁸O differencing line / Liu d18O
COL_BROWN_400 = '#8D6E63'   # KCM bars / Co secondary
COL_BROWN_200 = '#BCAAA4'   # δ¹⁸O fill

# ── Charcoal / blue-grey (axes, models, reference) ──
COL_BG_900    = '#263238'   # darkest — median lines, text
COL_BG_800    = '#37474F'   # PMIP4 bars, scatter points
COL_BG_700    = '#455A64'   # Co lines, d18O in c-panels
COL_BG_600    = '#546E7A'   # stat annotations
COL_BG_500    = '#607D8B'
COL_BG_400    = '#78909C'   # IQR fill (drip rate)
COL_BG_300    = '#90A4AE'   # reference lines, k_d
COL_BG_200    = '#B0BEC5'   # Co fill, event bands
COL_BG_100    = '#CFD8DC'
COL_BG_50     = '#ECEFF1'   # event shading

# ── Accent (used sparingly) ──
COL_RED_900   = '#B71C1C'   # trend lines, Tambora marker
COL_DORANGE   = '#D84315'   # TraCE bars, event markers
COL_GREEN_800 = '#2E7D32'   # baseflow marker
COL_INSOL     = '#90A4AE'   # blue-grey 300 — insolation curve

# ══════════════════════════════════════════════════════════════════════
# SEMANTIC ALIASES — use these in figure scripts
# ══════════════════════════════════════════════════════════════════════

# Fig 5 (drip rate hero)
COL_MEDIAN   = COL_BG_900
COL_IQR      = COL_BG_400
COL_D18O     = COL_BG_900
COL_D18O_C   = COL_BG_700     # d18O in c-panels
COL_LIU      = COL_BROWN_500  # Liu et al. high-res d18O
COL_AP       = COL_BG_300     # age-propagated envelope
COL_TREND    = COL_RED_900
COL_BASEFLOW = COL_GREEN_800
COL_MARKER   = COL_DORANGE

# Fig 6 (P reconstruction)
COL_HS4       = COL_TEAL_800
COL_HS4_IQR   = COL_TEAL_200
COL_HS4_90    = COL_TEAL_50
COL_D18O_P    = COL_BROWN_500
COL_D18O_FILL = COL_BROWN_200
COL_BAR_HS4   = COL_TEAL_600
COL_TRACE     = COL_DORANGE
COL_PMIP      = COL_BG_800
COL_KCM       = COL_BROWN_400

# Fig 3 (sensitivity)
COL_NI        = COL_TEAL_800
COL_NI_FILL   = COL_TEAL_200
COL_CO        = COL_BG_700
COL_CO_FILL   = COL_BG_200
COL_KD        = COL_BG_300

# Fig 4 (calibration)
COL_NI_LAB    = COL_TEAL_300
COL_CO_LAB    = COL_BG_400
COL_FIT       = COL_BG_900
COL_SCATTER   = COL_BG_800
COL_OBS       = COL_BG_900
COL_RECON     = COL_TEAL_600
COL_STAT      = COL_BG_600

# ══════════════════════════════════════════════════════════════════════
# EVENT BANDS — shared across Fig 5 and Fig 6
# ══════════════════════════════════════════════════════════════════════
EVENTS_FIG5 = {
    '8.2 ka':    {'c': 8200, 'hw': 250, 'col': COL_BG_200, 'a': 0.10},
    '5.2 ka':    {'c': 5200, 'hw': 200, 'col': COL_BROWN_200, 'a': 0.10},
    'Post-1750': {'c': 75,   'hw': 175, 'col': COL_BG_200, 'a': 0.13},
}

EVENTS_FIG6 = {
    '8.2 ka':    {'c': 8.2,   'hw': 0.25,  'col': COL_BG_50, 'a': 0.5},
    '5.2 ka':    {'c': 5.2,   'hw': 0.20,  'col': COL_BG_50, 'a': 0.5},
    'Post-1750': {'c': 0.075, 'hw': 0.175, 'col': COL_BG_50, 'a': 0.5},
}

# ══════════════════════════════════════════════════════════════════════
# HELPERS
# ══════════════════════════════════════════════════════════════════════

def style_ax(ax):
    """Apply NGeo spine and tick styling to an axes."""
    for sp in ax.spines.values():
        sp.set_linewidth(0.25)
    ax.tick_params(axis='both', which='major', width=0.25, length=2.5, direction='in')
    ax.tick_params(axis='both', which='minor', width=0.2,  length=1.2, direction='in')

def draw_events(axes, events):
    """Draw translucent event bands on a list of axes."""
    if not isinstance(axes, (list, tuple)):
        axes = [axes]
    for ev in events.values():
        c, hw = ev['c'], ev['hw']
        for ax in axes:
            ax.axvspan(c - hw, c + hw, color=ev['col'], alpha=ev['a'],
                       zorder=0.5, lw=0)

def make_heatmap_cmap(n=256):
    """Blue heatmap colourmap matched to Fig 2."""
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list('ngeo_fig2blue', [
        (1.0,  1.0,  1.0,  0.0),     # transparent white
        (0.97, 0.98, 1.00, 0.25),    # #F7FAFF — ice
        (0.84, 0.90, 0.96, 0.40),    # #D6E6F4 — pale
        (0.51, 0.73, 0.86, 0.60),    # #82BADB — sky
        (0.20, 0.51, 0.75, 0.80),    # #3282BE — mid
        (0.05, 0.34, 0.63, 0.92),    # #0C57A0 — dark
        (0.03, 0.19, 0.42, 1.0),     # #072F6B — navy
    ], N=n)

# ══════════════════════════════════════════════════════════════════════
# PAGE DIMENSIONS (mm → inches via /25.4)
# ══════════════════════════════════════════════════════════════════════
MM = 1 / 25.4   # multiply mm by this to get inches
SINGLE_COL = 89 * MM    # Nature single column
DOUBLE_COL = 180 * MM   # Nature double column
