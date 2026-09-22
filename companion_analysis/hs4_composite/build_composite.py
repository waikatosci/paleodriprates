#!/usr/bin/env python3
"""
build_composite.py -- stitches the eight hand-held DSLR photographs of the polished HS4 section into one
image on the pencil depth scale (Supplementary Figure 21), and cuts the axial slice used in Fig. 5c.

The photographs (C. Hu; CMYK JPEG, ~120 px/cm, not in this repository) are named by their depth range:
  B_top_35  A_35_70  D_70_97  C_98_131  E_128_157  F_157_191  G_189_222  H_208_242   (.jpg)
Set PHOTO_DIR to the folder holding them.

Steps
  1. CMYK -> sRGB through each file's embedded ICC profile.
  2. Find the pencil scale line in each frame (a robust fit to the darkest thin vertical feature) and its
     5 cm ticks (short strokes to the right of the line, darker than both the line's left side and the
     calcite beyond the tick). The first labelled tick of each frame is given in ANCHOR; the rest are
     numbered by a spacing search, and depth is fitted to image row (quadratic; residuals < 1 mm).
  3. Resample each frame to 100 px/cm with the scale line vertical at column XL, using the local
     px/cm along the line for both axes (square pixels) and the line's tilt for the perpendicular.
  4. Measure sub-mm offsets and colour ratios where frames of the same piece overlap (phase
     correlation; median ratios), apply them, and blend with weights that fall off towards frame edges.
     The black backdrop is masked (dark regions connected to the frame border); frames across a break
     are left at their own colour.
  5. Remove backdrop remnants at piece edges; cut the axial slice (-7 to -1 cm from the line).

Outputs (manuscript_figures/external/)
  HS4_composite.jpg / .json       100 px/cm, row 0 = 0 cm, scale line at column XL
  HS4_axial_slice.jpg / .json     6 cm wide slice, rows = depth, breaks white
"""
import glob
import json
import os
from pathlib import Path

import cv2
import numpy as np
from PIL import Image, ImageCms
from scipy.signal import find_peaks

REPO = Path(__file__).resolve().parents[2]
EXT = Path(os.environ.get('OUT_DIR', REPO / 'manuscript_figures' / 'external'))
PHOTO_DIR = Path(os.environ.get('PHOTO_DIR', '.'))
WORK = Path(os.environ.get('WORK_DIR', PHOTO_DIR))
ORDER = 'BADCEFGH'                                  # top to base
ANCHOR = {'B': (750, 10), 'A': (21, 35), 'D': (93, 70), 'C': (265, 100),    # (row of first labelled tick,
          'E': (334, 130), 'F': (434, 160), 'G': (14, 190), 'H': (185, 210)}  #  its depth in cm)
PX, XL, W_OUT, DEPTH_MAX = 100, 1650, 2450, 246
SLICE = (-7.0, -1.0)
np.random.seed(20260923)


def to_srgb(c):
    src = glob.glob(str(PHOTO_DIR / f'{c}_*.jpg'))[0]
    im = Image.open(src)
    if im.mode == 'CMYK':
        prof = ImageCms.ImageCmsProfile(__import__('io').BytesIO(im.info['icc_profile']))
        im = ImageCms.profileToProfile(im, prof, ImageCms.createProfile('sRGB'), outputMode='RGB')
    return cv2.cvtColor(np.asarray(im.convert('RGB')), cv2.COLOR_RGB2BGR)


def line_and_ticks(im):
    g = cv2.cvtColor(im, cv2.COLOR_BGR2GRAY).astype(np.float32)
    h, w = g.shape
    hp = cv2.blur(g, (41, 1)) - g
    hp[g < 60] = 0
    vv = cv2.blur(hp, (1, 301))
    ys = np.arange(150, h - 150, 50)
    xs = np.array([np.argmax(vv[y, w // 3:]) + w // 3 for y in ys], float)
    ok = np.array([(g[y - 100:y + 100, int(x)] > 60).mean() > 0.9 for y, x in zip(ys, xs)])
    best = None
    for _ in range(300):
        a, b = np.random.choice(np.where(ok)[0], 2, replace=False)
        p = np.polyfit(ys[[a, b]], xs[[a, b]], 1)
        inl = ok & (np.abs(xs - np.polyval(p, ys)) < 4)
        if best is None or inl.sum() > best[1].sum():
            best = (p, inl)
    line = np.polyfit(ys[best[1]], xs[best[1]], 2)
    xl = np.round(np.polyval(line, np.arange(h))).astype(int)
    tick = np.array([g[y, xl[y] + 5:xl[y] + 28].mean() for y in range(h)])
    left = np.array([g[y, xl[y] - 40:xl[y] - 8].mean() for y in range(h)])
    far = np.array([g[y, xl[y] + 45:xl[y] + 75].mean() for y in range(h)])
    d = np.convolve(np.minimum(left, far) - tick, np.ones(5) / 5, 'same')
    pk1, _ = find_peaks(d, height=8, distance=250)
    # second detector: tick darkness against a vertically smoothed reference
    ref = cv2.blur(g, (1, 61))
    prof = np.array([-(g[y, xl[y] + 6:xl[y] + 30].mean()) for y in range(h)])
    refp = np.array([-(ref[y, xl[y] + 6:xl[y] + 30].mean()) for y in range(h)])
    pk2, _ = find_peaks(prof - refp, height=6, distance=400, prominence=6)
    return line, np.r_[pk1, pk2], np.r_[d[pk1], np.full(len(pk2), 30.0)]


def depth_fit(peaks, hts, anchor):
    ya, da = anchor
    best = None
    for s in np.arange(540, 620, 0.5):
        k = np.round((peaks - ya) / s)
        r = peaks - (ya + s * k)
        m = (np.abs(r) < 20) & (hts > 20)
        sc = len(set(k[m])) * 100 - np.abs(r[m]).mean()
        if best is None or sc > best[0]:
            best = (sc, s)
    s = best[1]
    k = np.round((peaks - ya) / s)
    m = (np.abs(peaks - (ya + s * k)) < 20) & (hts > 20)
    ks = sorted(set(k[m]))
    Y = np.array([peaks[m][k[m] == q].mean() for q in ks])
    D = da + 5 * np.array(ks)
    fit = np.polyfit(Y, D, 2 if len(Y) >= 5 else 1)
    return fit, np.abs(D - np.polyval(fit, Y)).max()


def remap(im, line, fit):
    h = im.shape[0]
    yy = np.arange(h, dtype=float)
    dd = np.polyval(fit, yy)
    r0, r1 = int(np.ceil(dd[0] * PX)), int(np.floor(dd[-1] * PX))
    d = np.arange(r0, r1) / PX
    ys = np.interp(d, dd, yy)
    k = 1 / np.polyval(np.polyder(fit), ys)
    xls = np.polyval(line, ys)
    sl = np.polyval(np.polyder(line), ys)
    nrm = np.sqrt(1 + sl ** 2)
    u = (np.arange(W_OUT) - XL) / PX
    mx = (xls[:, None] + u[None, :] * (k / nrm)[:, None]).astype(np.float32)
    my = (ys[:, None] + u[None, :] * (k * -sl / nrm)[:, None]).astype(np.float32)
    return cv2.remap(im, mx, my, cv2.INTER_LINEAR, borderMode=cv2.BORDER_CONSTANT, borderValue=0), r0


def stone_mask(im):
    gray = cv2.cvtColor(im, cv2.COLOR_BGR2GRAY)
    dark = (gray < 65).astype(np.uint8)
    bg = cv2.morphologyEx(dark, cv2.MORPH_OPEN, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (35, 35)))
    _, lab = cv2.connectedComponents(bg)
    edge = set(np.unique(np.concatenate([lab[0], lab[-1], lab[:, 0], lab[:, -1]]))) - {0}
    bg = cv2.dilate(np.isin(lab, list(edge)).astype(np.uint8), np.ones((9, 9)))
    inb = cv2.erode((im.max(2) > 0).astype(np.uint8), np.ones((15, 15)))
    return ((1 - bg) * inb).astype(np.float32)


def overlap_offsets(a, ra, b, rb, ma, mb):
    r0, r1 = max(ra, rb), min(ra + len(a), rb + len(b))
    if r1 - r0 < 50:
        return None
    A, B = a[r0 - ra:r1 - ra], b[r0 - rb:r1 - rb]
    m = (ma[r0 - ra:r1 - ra] > 0) & (mb[r0 - rb:r1 - rb] > 0)
    if m.sum() < 1e5:
        return None
    ratio = np.array([np.median(A[..., i][m].astype(float)) / np.median(B[..., i][m].astype(float)) for i in range(3)])
    # offset of B relative to A: normalised cross-correlation of the high-passed overlap, searched over
    # +-15 px (1.5 mm), with a parabolic sub-pixel peak
    frac = m.mean(0)
    cols = np.where(frac > min(0.9, 0.9 * frac.max()))[0]
    x0, x1 = cols.min() + 20, cols.max() - 20
    ga = cv2.cvtColor(A, cv2.COLOR_BGR2GRAY).astype(np.float32)
    gb = cv2.cvtColor(B, cv2.COLOR_BGR2GRAY).astype(np.float32)
    ga -= cv2.GaussianBlur(ga, (0, 0), 20)
    gb -= cv2.GaussianBlur(gb, (0, 0), 20)
    S = 15
    tpl = gb[S:-S, x0:x1]
    res = cv2.matchTemplate(ga[:, x0 - S:x1 + S], tpl, cv2.TM_CCOEFF_NORMED)
    iy, ix = np.unravel_index(np.argmax(res), res.shape)
    def sub(v, i):
        if 0 < i < len(v) - 1:
            den = v[i - 1] - 2 * v[i] + v[i + 1]
            return i + (0.5 * (v[i - 1] - v[i + 1]) / den if den else 0)
        return float(i)
    dy = S - sub(res[:, ix], iy)          # B content sits at A + (dx, dy)
    dx = S - sub(res[iy, :], ix)
    return (dx, dy), ratio, res.max(), res[S, S]


def main():
    frames = {}
    for c in ORDER:
        im = to_srgb(c)
        line, pk, ht = line_and_ticks(im)
        fit, res = depth_fit(pk, ht, ANCHOR[c])
        rm, r0 = remap(im, line, fit)
        frames[c] = dict(im=rm, r0=r0, mask=stone_mask(rm))
        print(f'{c}: {np.polyval(fit, 0):.1f}-{np.polyval(fit, im.shape[0] - 1):.1f} cm, max tick residual {res * 10:.2f} mm')
    # offsets and colour ratios between consecutive frames of the same piece
    gain = {ORDER[0]: np.ones(3)}
    for a, b in zip(ORDER, ORDER[1:]):
        o = overlap_offsets(frames[a]['im'], frames[a]['r0'], frames[b]['im'], frames[b]['r0'],
                            frames[a]['mask'], frames[b]['mask'])
        if o is None:                      # frames across a break: own colour, no shift
            gain[b] = np.ones(3)
            print(f'{a}/{b}: across a break')
            continue
        (dx, dy), ratio, ncc, ncc0 = o
        gain[b] = gain[a] * ratio
        t = np.float32([[1, 0, -dx], [0, 1, -dy]])
        fb = frames[b]
        fb['im'] = cv2.warpAffine(fb['im'], t, (W_OUT, len(fb['im'])), borderValue=0)
        fb['mask'] = cv2.warpAffine(fb['mask'], t, (W_OUT, len(fb['mask'])), borderValue=0)
        print(f'{a}/{b}: offset {dx / PX * 10:.2f}, {dy / PX * 10:.2f} mm (match {ncc0:.2f} -> {ncc:.2f}); colour ratio {np.round(ratio, 3)}')
    H = DEPTH_MAX * PX
    acc = np.zeros((H, W_OUT, 3), np.float32)
    wsum = np.zeros((H, W_OUT), np.float32)
    for c in ORDER:
        f = frames[c]
        n = len(f['im'])
        ramp = np.minimum(np.arange(n), np.arange(n)[::-1]).astype(np.float32) + 1
        wt = f['mask'] * np.minimum(ramp, 150)[:, None]
        wt = cv2.GaussianBlur(wt, (0, 0), 3) * f['mask']
        acc[f['r0']:f['r0'] + n] += f['im'].astype(np.float32) * gain[c][None, None, :] * wt[..., None]
        wsum[f['r0']:f['r0'] + n] += wt
    out = np.full((H, W_OUT, 3), 255, np.float32)
    k = wsum > 1e-2
    out[k] = acc[k] / wsum[k][:, None]
    out = np.clip(out, 0, 255).astype(np.uint8)
    # backdrop remnants at piece edges: near-black components touching the white background
    g = cv2.cvtColor(out, cv2.COLOR_BGR2GRAY)
    white = out.min(2) >= 250
    n, lab, st, _ = cv2.connectedComponentsWithStats((g < 70).astype(np.uint8), 8)
    touch = np.unique(lab[cv2.dilate(white.astype(np.uint8), np.ones((5, 5))) > 0])
    kill = np.zeros(n, bool)
    for i in touch[touch > 0]:
        m = lab == i
        if st[i, 4] > 300 and g[m].mean() < 45:
            kill[i] = True
    out[kill[lab]] = 255
    cv2.imwrite(str(EXT / 'HS4_composite.jpg'), out, [cv2.IMWRITE_JPEG_QUALITY, 90])
    json.dump({'px_per_cm': PX, 'row0_depth_cm': 0.0, 'depth_line_column_px': XL, 'background_rgb': [255, 255, 255]},
              open(EXT / 'HS4_composite.json', 'w'))
    s = out[:, XL + int(SLICE[0] * PX):XL + int(SLICE[1] * PX)].copy()
    s[s.min(2) >= 250] = 255
    cv2.imwrite(str(EXT / 'HS4_axial_slice.jpg'), s, [cv2.IMWRITE_JPEG_QUALITY, 92])
    json.dump({'px_per_cm': PX, 'row0_depth_cm': 0.0, 'x_from_line_cm': list(SLICE), 'break_fill_rgb': [255, 255, 255]},
              open(EXT / 'HS4_axial_slice.json', 'w'))
    print('saved HS4_composite and HS4_axial_slice')


if __name__ == '__main__':
    main()
