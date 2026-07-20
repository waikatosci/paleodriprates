#!/usr/bin/env python3
"""
generate_figures.py — regenerate every main-text figure for NCOMMS-26-041445-T
from the single canonical source, HS4_SourceData.xlsx.

    python3 generate_figures.py                 # all figures
    python3 generate_figures.py 3 5             # only Figs 3 and 5
    python3 generate_figures.py --check         # preflight only, render nothing
    python3 generate_figures.py --list          # what exists and what is blocked

WHAT THIS IS
    Every figure script here reads its data from HS4_SourceData.xlsx and nothing else,
    except for a short list of bulk arrays that do not belong in a spreadsheet (a 131 MB
    basemap, a 248 MB Monte-Carlo ensemble, a 119k-cell PDF grid). Those are listed in
    external/README.md with checksums.

    Every plotted value has a single home in the workbook, so the same measurement
    cannot exist in two places and drift between them. Everything below comes

PREFLIGHT
    --check reports missing inputs and fonts. Render finals on one machine with Arial
    installed, then check with `pdffonts` that exactly one family comes back.

EXIT CODES
    0 all requested figures rendered
    1 something failed or is blocked
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'HS4_SourceData.xlsx'

# figure id -> (script, human name, extra files it needs beyond the workbook)
FIGURES = {
    '2': ('figures/Fig2_site.py', 'Site map + monitoring',
          ['external/precip.mon.mean.nc', 'external/SR_LR.tif']),
    '3': ('figures/Fig3_sensitivity.py', 'Forward-model sensitivity', []),
    '4': ('figures/Fig4_calibration.py', 'Calibration', []),
    '5': ('figures/Fig5_record.py', 'Holocene record',
          ['external/pdf_heatmap_hr.json', 'external/drip_rate_summary_hr.csv']),
    '6': ('figures/Fig6_P_comparison.py', 'Precipitation comparison', []),
    '7': ('figures/Fig7_events.py', 'Event detection',
          ['external/pdf_heatmap_hr.json', 'external/drip_rate_summary_hr.csv',
           'external/drip_rate_summary_lr.csv', 'external/HS4_Zhu2017_IRMsoft_flux.csv']),
}

# Fig 1 is a conceptual diagram with no data; it is not generated here.
NOT_GENERATED = {'1': 'conceptual diagram — no data, not generated from source'}


def _c(s, colour):
    if not sys.stdout.isatty():
        return s
    return {'r': f'\033[31m{s}\033[0m', 'g': f'\033[32m{s}\033[0m',
            'y': f'\033[33m{s}\033[0m', 'b': f'\033[1m{s}\033[0m'}.get(colour, s)


def check_fonts():
    """Report which sans family will actually be used. See the module docstring."""
    try:
        from matplotlib import font_manager
    except ImportError:
        return None, ['matplotlib not importable']
    have = {f.name for f in font_manager.fontManager.ttflist}
    stack = ['Helvetica', 'Arial', 'Liberation Sans', 'DejaVu Sans']
    resolved = next((n for n in stack if n in have), None)
    notes = []
    if resolved not in ('Helvetica', 'Arial'):
        notes.append(f'Arial/Helvetica NOT installed — figures will render in '
                     f'{resolved!r}. Fine for a working render; NOT for finals.')
    if resolved == 'DejaVu Sans':
        notes.append('DejaVu is NOT metric-compatible with Arial: label metrics and '
                     'therefore layout will differ, not just glyph shapes.')
    try:
        import matplotlib
        if matplotlib.rcParams.get('mathtext.fontset') == 'dejavusans':
            notes.append("mathtext.fontset is matplotlib's 'dejavusans' default: maths "
                         "renders in DejaVu even when Arial is present. "
                         "check ngeo_style.py font settings.")
    except Exception:
        pass
    return resolved, notes


def preflight(ids):
    ok = True
    print(_c('Preflight', 'b'))
    if SOURCE.exists():
        print(f'  {_c("OK", "g"):>6}  HS4_SourceData.xlsx  ({SOURCE.stat().st_size/1e6:.1f} MB)')
    else:
        print(f'  {_c("MISSING", "r"):>6}  HS4_SourceData.xlsx  <- nothing can run without this')
        ok = False
    if not (ROOT / 'ngeo_style.py').exists():
        print(f'  {_c("MISSING", "r"):>6}  ngeo_style.py')
        ok = False

    for mod in ('pandas', 'numpy', 'matplotlib', 'scipy', 'openpyxl'):
        try:
            __import__(mod)
        except ImportError:
            print(f'  {_c("MISSING", "r"):>6}  python package: {mod}')
            ok = False
    for fid in ids:
        script, name, extras = FIGURES[fid]
        if not (ROOT / script).exists():
            print(f'  {_c("MISSING", "r"):>6}  {script}')
            ok = False
        for e in extras:
            if not (ROOT / e).exists():
                print(f'  {_c("MISSING", "r"):>6}  {e}   (needed by Fig {fid})')
                ok = False

    resolved, notes = check_fonts()
    if resolved:
        tag = 'OK' if resolved in ('Helvetica', 'Arial') else 'WARN'
        print(f'  {_c(tag, "g" if tag == "OK" else "y"):>6}  fonts: will render in {resolved!r}')
    for n in notes:
        print(f'          {_c("!", "y")} {n}')
    return ok


def run_one(fid):
    script, name, _ = FIGURES[fid]
    print(f'\n{_c(f"Fig {fid}", "b")} — {name}   [{script}]')
    r = subprocess.run([sys.executable, str(ROOT / script)], cwd=ROOT,
                       capture_output=True, text=True,
                       env=dict(os.environ, MPLBACKEND='Agg'))
    for line in r.stdout.splitlines():
        if line.strip():
            print(f'    {line}')
    if r.returncode != 0:
        print(_c(f'    FAILED (exit {r.returncode})', 'r'))
        tail = [l for l in r.stderr.splitlines() if l.strip()][-6:]
        for line in tail:
            print(f'    {_c(line, "r")}')
        return False
    print(_c('    OK', 'g'))
    return True


def main():
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument('ids', nargs='*', help='figure numbers (default: all)')
    ap.add_argument('--check', action='store_true', help='preflight only')
    ap.add_argument('--list', action='store_true', help='list figures and status')
    a = ap.parse_args()

    if a.list:
        print(_c('Figures', 'b'))
        for fid, (script, name, extras) in FIGURES.items():
            miss = [e for e in extras if not (ROOT / e).exists()]
            have_script = (ROOT / script).exists()
            if not have_script:
                st = _c('BLOCKED', 'r') + '  script not present'
            elif miss:
                st = _c('BLOCKED', 'r') + '  missing: ' + ', '.join(os.path.basename(m) for m in miss)
            else:
                st = _c('ready', 'g')
            print(f'  Fig {fid}  {name:32s} {st}')
        for fid, why in NOT_GENERATED.items():
            print(f'  Fig {fid}  {"":32s} {_c("n/a", "y")}  {why}')
        return 0

    ids = a.ids or [f for f in FIGURES if (ROOT / FIGURES[f][0]).exists()]
    bad = [i for i in ids if i not in FIGURES]
    if bad:
        print(f'unknown figure(s): {bad}. Known: {sorted(FIGURES)}')
        return 1

    ok = preflight(ids)
    if a.check:
        return 0 if ok else 1
    if not ok:
        print(_c('\nPreflight failed — fix the above, or run with a subset.', 'r'))
        return 1

    results = {}
    for fid in ids:
        results[fid] = run_one(fid)

    print('\n' + '─' * 60)
    for fid, good in results.items():
        print(f'  Fig {fid}  {_c("OK", "g") if good else _c("FAILED", "r")}')
    n_ok = sum(results.values())
    print(f'  {n_ok}/{len(results)} figures')
    if n_ok == len(results):
        print(_c('\n  Reminder: check fonts before submitting finals —', 'y'))
        print('    pdffonts <figure>.pdf     # expect exactly ONE family')
    return 0 if n_ok == len(results) else 1


if __name__ == '__main__':
    sys.exit(main())
