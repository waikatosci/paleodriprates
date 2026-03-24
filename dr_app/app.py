"""
Drip Rate Estimator - Flask Web Application
============================================
Provides a browser-based interface for the speleothem drip rate
reconstruction pipeline (replaces the Excel-based workflow).
"""

import os, sys, json, threading, time, traceback, pickle, copy, datetime
from flask import Flask, request, jsonify, send_file, Response
import numpy as np
import pandas as pd

app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 50 * 1024 * 1024  # 50 MB upload limit
UPLOAD_FOLDER = os.path.join(os.path.dirname(__file__), 'uploads')
OUTPUT_FOLDER = os.path.join(os.path.dirname(__file__), 'outputs')
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# ── Global run state ────────────────────────────────────────────────────────
run_state = {
    'running': False,
    'progress': 0,
    'stage': '',
    'log': [],
    'error': None,
    'done': False,
    'outputs': [],
}

def log(msg):
    run_state['log'].append(msg)
    print(msg)


# ── Routes ───────────────────────────────────────────────────────────────────

HTML = r'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<meta http-equiv="Cache-Control" content="no-cache, no-store, must-revalidate">
<meta http-equiv="Pragma" content="no-cache">
<meta http-equiv="Expires" content="0">
<title>Drip Rate Estimator (v44)</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/4.4.1/chart.umd.min.js"></script>
<style>
  :root {
    --bg:       #0d1117;
    --surface:  #161b22;
    --border:   #2a3441;
    --accent:   #4cc9a0;
    --accent2:  #f7a440;
    --muted:    #8fa4b5;
    --text:     #dce6ef;
    --texthi:   #f0f4f8;
    --danger:   #e05c5c;
    --radius:   10px;
    --mono:     Arial, Helvetica, sans-serif;
    --serif:    Arial, Helvetica, sans-serif;
  }

  * { box-sizing: border-box; margin: 0; padding: 0; }

  body {
    background: var(--bg);
    color: var(--text);
    font-family: var(--mono);
    font-size: 14px;
    min-height: 100vh;
  }

  /* ── Layout ── */
  .shell {
    display: grid;
    grid-template-columns: 280px 1fr;
    grid-template-rows: 56px 1fr auto;
    height: 100vh;
    overflow: hidden;
  }

  header {
    grid-column: 1 / -1;
    background: var(--surface);
    border-bottom: 1px solid var(--border);
    display: flex;
    align-items: center;
    padding: 0 24px;
    gap: 12px;
  }

  header .logo {
    font-family: var(--serif);
    font-size: 20px;
    font-weight: 600;
    color: var(--accent);
    letter-spacing: -0.5px;
  }

  header .subtitle {
    color: var(--muted);
    font-size: 11px;
    border-left: 1px solid var(--border);
    padding-left: 12px;
    margin-left: 4px;
  }

  nav {
    background: var(--surface);
    border-right: 1px solid var(--border);
    overflow-y: auto;
    padding: 8px 0;
  }

  .nav-section {
    padding: 16px 16px 4px;
    font-size: 10px;
    text-transform: uppercase;
    letter-spacing: 1.5px;
    color: var(--muted);
  }

  .nav-item {
    display: flex;
    align-items: center;
    gap: 10px;
    padding: 9px 16px;
    cursor: pointer;
    color: var(--text);
    border-left: 3px solid transparent;
    transition: all 0.15s;
    user-select: none;
  }
  .nav-item:hover { background: rgba(76,201,160,0.06); }
  .nav-item.active {
    color: var(--accent);
    border-left-color: var(--accent);
    background: rgba(76,201,160,0.08);
  }
  .nav-item .icon { font-size: 15px; width: 18px; text-align: center; }
  .nav-item .badge {
    margin-left: auto;
    background: var(--accent);
    color: var(--bg);
    border-radius: 10px;
    padding: 1px 7px;
    font-size: 10px;
  }
  .nav-item .badge.warn { background: var(--accent2); }

  main {
    overflow-y: auto;
    padding: 28px 32px;
  }

  /* ── Panels ── */
  .panel { display: none; }
  .panel.active { display: block; }

  .page-title {
    font-family: var(--serif);
    font-size: 24px;
    font-weight: 300;
    color: var(--texthi);
    margin-bottom: 6px;
  }
  .page-desc {
    color: var(--muted);
    margin-bottom: 24px;
    line-height: 1.6;
    font-size: 12px;
  }

  /* ── Cards ── */
  .card {
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: var(--radius);
    padding: 20px 24px;
    margin-bottom: 16px;
  }
  .card-title {
    font-size: 12px;
    text-transform: uppercase;
    letter-spacing: 1.2px;
    color: var(--accent);
    margin-bottom: 16px;
    display: flex;
    align-items: center;
    gap: 8px;
  }

  /* ── Forms ── */
  .form-grid {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(200px, 1fr));
    gap: 14px;
  }
  .form-grid.wide { grid-template-columns: repeat(auto-fill, minmax(280px, 1fr)); }

  .field { display: flex; flex-direction: column; gap: 5px; }
  .field label { font-size: 12px; color: var(--muted); }
  .field input, .field select {
    background: var(--bg);
    border: 1px solid var(--border);
    border-radius: 6px;
    color: var(--texthi);
    font-family: var(--mono);
    font-size: 13px;
    padding: 8px 10px;
    transition: border-color 0.15s;
    width: 100%;
  }
  .field input:focus, .field select:focus {
    outline: none;
    border-color: var(--accent);
  }
  .field select option { background: var(--surface); }

  /* ── Upload zones ── */
  .upload-grid {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(220px, 1fr));
    gap: 12px;
    margin-bottom: 8px;
  }

  .dropzone {
    border: 2px dashed var(--border);
    border-radius: var(--radius);
    padding: 20px 14px;
    text-align: center;
    cursor: pointer;
    transition: all 0.2s;
    position: relative;
  }
  .dropzone:hover, .dropzone.drag { border-color: var(--accent); background: rgba(76,201,160,0.04); }
  .dropzone.ready { border-color: var(--accent); border-style: solid; }
  .dropzone .dz-icon { font-size: 22px; margin-bottom: 6px; }
  .dropzone .dz-label { font-size: 12px; color: var(--muted); margin-bottom: 4px; }
  .dropzone .dz-name { font-size: 12px; color: var(--accent); display: none; }
  .dropzone.ready .dz-name { display: block; }
  .dropzone.ready .dz-hint { display: none; }
  .dropzone .dz-hint { font-size: 10px; color: var(--muted); }
  .dropzone input[type=file] {
    position: absolute; inset: 0; opacity: 0; cursor: pointer; width: 100%; height: 100%;
  }

  /* column selector that appears after upload */
  .col-selector { margin-top: 10px; display: none; }
  .col-selector.show { display: block; }
  .col-selector label { font-size: 11px; color: var(--muted); display: block; margin-bottom: 3px; }
  .col-selector select {
    width: 100%;
    background: var(--bg);
    border: 1px solid var(--border);
    border-radius: 6px;
    color: var(--texthi);
    font-family: var(--mono);
    font-size: 11px;
    padding: 6px 8px;
    margin-bottom: 6px;
  }

  /* ── Buttons ── */
  .btn {
    display: inline-flex;
    align-items: center;
    gap: 7px;
    padding: 10px 20px;
    border-radius: 7px;
    border: none;
    font-family: var(--mono);
    font-size: 12px;
    cursor: pointer;
    transition: all 0.15s;
    font-weight: 500;
  }
  .btn-primary {
    background: var(--accent);
    color: var(--bg);
  }
  .btn-primary:hover { filter: brightness(1.1); }
  .btn-primary:disabled { opacity: 0.4; cursor: not-allowed; }
  .btn-ghost {
    background: transparent;
    border: 1px solid var(--border);
    color: var(--text);
  }
  .btn-ghost:hover { border-color: var(--accent); color: var(--accent); }
  .btn-danger {
    background: transparent;
    border: 1px solid var(--danger);
    color: var(--danger);
  }

  /* ── Mode toggle buttons (manual/CSV) ── */
  .mode-btn {
    background: var(--surface);
    border: 1px solid var(--border);
    color: var(--muted);
    border-radius: 6px;
    padding: 4px 12px;
    font-family: var(--mono);
    font-size: 12px;
    cursor: pointer;
    transition: background 0.15s, color 0.15s, border-color 0.15s;
  }
  .mode-btn.active {
    background: var(--accent);
    color: #0d1117;
    border-color: var(--accent);
    font-weight: 500;
  }
  .mode-btn:hover:not(.active) {
    border-color: var(--accent);
    color: var(--accent);
  }

  /* ── Progress ── */
  .progress-wrap {
    background: var(--bg);
    border-radius: 99px;
    height: 6px;
    overflow: hidden;
    margin: 12px 0 6px;
  }
  .progress-bar {
    height: 100%;
    background: linear-gradient(90deg, var(--accent), #7be0c0);
    border-radius: 99px;
    transition: width 0.4s ease;
    width: 0%;
  }
  .progress-label {
    display: flex;
    justify-content: space-between;
    font-size: 11px;
    color: var(--muted);
  }

  /* ── Log ── */
  .log-box {
    background: var(--bg);
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 14px;
    height: 220px;
    overflow-y: auto;
    font-size: 11px;
    line-height: 1.7;
    color: var(--muted);
  }
  .log-box .log-line { }
  .log-box .log-ok  { color: var(--accent); }
  .log-box .log-err { color: var(--danger); }

  /* ── Status pill ── */
  .status-pill {
    display: inline-flex;
    align-items: center;
    gap: 6px;
    padding: 4px 12px;
    border-radius: 99px;
    font-size: 11px;
    font-weight: 500;
    margin-left: auto;
  }
  .status-pill.idle    { background: rgba(110,127,141,0.15); color: var(--muted); }
  .status-pill.running { background: rgba(76,201,160,0.15);  color: var(--accent); }
  .status-pill.done    { background: rgba(76,201,160,0.2);   color: var(--accent); }
  .status-pill.error   { background: rgba(224,92,92,0.15);   color: var(--danger); }
  .dot { width: 7px; height: 7px; border-radius: 50%; background: currentColor; }
  .dot.pulse { animation: pulse 1.2s ease infinite; }
  @keyframes pulse { 0%,100% { opacity:1; } 50% { opacity:0.3; } }

  /* ── Chart ── */
  .chart-wrap {
    background: var(--bg);
    border-radius: 8px;
    padding: 16px;
    position: relative;
    height: 340px;
  }

  /* ── Output file list ── */
  .output-list { display: flex; flex-direction: column; gap: 8px; }
  .output-item {
    display: flex;
    align-items: center;
    gap: 12px;
    padding: 12px 16px;
    background: var(--bg);
    border: 1px solid var(--border);
    border-radius: 8px;
  }
  .output-item .oi-icon { font-size: 18px; }
  .output-item .oi-name { flex: 1; font-size: 12px; color: var(--texthi); }
  .output-item .oi-desc { font-size: 11px; color: var(--muted); }

  /* ── Info callout ── */
  .callout {
    background: rgba(76,201,160,0.06);
    border: 1px solid rgba(76,201,160,0.2);
    border-radius: 8px;
    padding: 12px 16px;
    font-size: 12px;
    color: var(--text);
    line-height: 1.6;
    margin-bottom: 16px;
  }
  .callout.warn {
    background: rgba(247,164,64,0.06);
    border-color: rgba(247,164,64,0.25);
    color: var(--accent2);
  }

  hr.divider { border: none; border-top: 1px solid var(--border); margin: 20px 0; }

  /* ── Scrollbar ── */
  ::-webkit-scrollbar { width: 6px; }
  ::-webkit-scrollbar-track { background: transparent; }
  ::-webkit-scrollbar-thumb { background: var(--border); border-radius: 3px; }

  /* ── Tooltips ── */
  .has-tip { position: relative; display: flex; flex-direction: column; gap: 5px; }
  .tip-icon {
    display: inline-flex; align-items: center; justify-content: center;
    width: 15px; height: 15px; border-radius: 50%;
    background: var(--border); color: var(--muted);
    font-size: 9px; font-style: normal; cursor: help;
    flex-shrink: 0; margin-left: 5px; vertical-align: middle;
    transition: background 0.15s;
  }
  .tip-icon:hover { background: var(--accent); color: var(--bg); }
  .tip-label { display: flex; align-items: center; }
  .tooltip-box {
    display: none;
    position: absolute;
    top: calc(100% + 8px);
    left: 50%;
    transform: translateX(-50%);
    right: auto;
    max-width: min(640px, calc(100vw - 40px));
    width: 640px;
    background: #1c2635;
    overflow-y: auto;
    max-height: 80vh;
    border: 1px solid var(--accent);
    border-radius: 8px;
    padding: 12px 14px;
    font-size: 12px;
    line-height: 1.7;
    color: var(--text);
    z-index: 999;
    box-shadow: 0 8px 24px rgba(0,0,0,0.5);
    pointer-events: none;
  }
  .tooltip-box .tb-title {
    font-size: 11px; color: var(--accent); font-weight: 500;
    margin-bottom: 6px; text-transform: uppercase; letter-spacing: 0.8px;
  }
  .tooltip-box code {
    display: block; background: var(--bg); border-radius: 4px;
    padding: 6px 8px; margin: 6px 0; font-size: 11px;
    color: var(--accent2); white-space: pre-wrap;
  }
  .has-tip:hover .tooltip-box { display: block; }


  /* ── Footer ── */
  .app-footer {
    grid-column: 1 / -1;
    background: var(--surface);
    border-top: 1px solid var(--border);
    padding: 10px 24px;
    font-size: 11px;
    color: var(--muted);
    line-height: 1.7;
    display: flex;
    align-items: flex-start;
    gap: 10px;
    flex-wrap: wrap;
  }
  .footer-badge {
    background: rgba(247,164,64,0.15);
    color: var(--accent2);
    border: 1px solid rgba(247,164,64,0.3);
    border-radius: 4px;
    padding: 2px 8px;
    font-size: 10px;
    white-space: nowrap;
    flex-shrink: 0;
    align-self: center;
  }
  .footer-text { flex: 1; min-width: 280px; }
  .footer-text em  { color: var(--text); }
  .footer-text strong { color: var(--texthi); }
  .app-footer a { color: var(--accent); text-decoration: none; }
  .app-footer a:hover { text-decoration: underline; }


  /* ── Sidebar image ── */
  .nav-image-wrap {
    position: relative;
    width: 100%;
    overflow: hidden;
    border-bottom: 1px solid var(--border);
  }
  .nav-image-wrap img {
    width: 100%;
    display: block;
    object-fit: cover;
    filter: brightness(0.82) saturate(0.9);
    transition: filter 0.3s;
  }
  .nav-image-wrap:hover img { filter: brightness(0.95) saturate(1.1); }

  /* ── Drip animation (active during model run) ── */
  .nav-image-wrap.running img {
    animation: dripPulse 2.2s ease-in-out infinite;
  }
  @keyframes dripPulse {
    0%   { filter: brightness(0.75) saturate(0.9); }
    50%  { filter: brightness(1.1)  saturate(1.3); }
    100% { filter: brightness(0.75) saturate(0.9); }
  }
  .nav-image-wrap.running .nav-image-caption {
    border-bottom: 2px solid var(--accent);
    transition: border-color 0.3s;
  }

  .nav-image-caption {
    position: absolute;
    bottom: 0; left: 0; right: 0;
    background: linear-gradient(transparent, rgba(0,0,0,0.82));
    padding: 18px 10px 8px;
    font-size: 9px;
    color: rgba(255,255,255,0.65);
    line-height: 1.5;
  }
  .nav-image-caption a {
    color: var(--accent);
    text-decoration: none;
  }
  .nav-image-caption a:hover { text-decoration: underline; }

</style>
</head>
<body>
<div class="shell">

  <!-- ── Header ── -->
  <header>
    <span class="logo">⬡ DripRate</span>
    <span class="subtitle">Speleothem Palaeoclimate Reconstruction</span>
    <span id="statusPill" class="status-pill idle">
      <span class="dot"></span> Idle
    </span>
  </header>

  <!-- ── Sidebar ── -->
  <nav>
    <!-- Sidebar splash image -->
    <div class="nav-image-wrap">
      <img src="data:image/jpeg;base64,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" alt="Cave drip splash on stalagmite tip">
      <div class="nav-image-caption">
        © Garry Smith<br>
        <a href="https://www.geochemicalperspectivesletters.org/article1824/"
           target="_blank" title="Hartland & Zitoun, Geochemical Perspectives Letters">
          Hartland &amp; Zitoun, Geochem. Persp. Lett.
        </a>
      </div>
    </div>

    <div class="nav-section">Workflow</div>
    <div class="nav-item active" onclick="showPanel('data')">
      <span class="icon">📂</span> Data Input
      <span id="badgeData" class="badge warn" style="display:none">!</span>
    </div>
    <div class="nav-item" onclick="showPanel('params')">
      <span class="icon">⚙️</span> Model Parameters
    </div>
    <div class="nav-item" onclick="showPanel('output')">
      <span class="icon">🎛️</span> Output Options
    </div>
    <div class="nav-item" onclick="showPanel('run')">
      <span class="icon">▶</span> Run Model
    </div>
    <div class="nav-section">Info</div>
    <div class="nav-item" onclick="showPanel('about')">
      <span class="icon">ℹ️</span> About
    </div>
    <div class="nav-section">Results</div>
    <div class="nav-item" onclick="showPanel('results')">
      <span class="icon">📈</span> Results
      <span id="badgeResults" class="badge" style="display:none">New</span>
    </div>
    <div class="nav-item" onclick="showPanel('downloads')">
      <span class="icon">⬇</span> Downloads
    </div>
  </nav>

  <!-- ── Main content ── -->
  <main>

    <!-- DATA INPUT -->
    <div id="panel-data" class="panel active">
      <div class="page-title">Data Input</div>
      <p class="page-desc">Upload your CSV files and map columns to the expected variables. All files should have a header row.</p>

      <div class="callout">
        Expected format: comma-separated CSV with a header row. Depth in <strong>cm</strong>, age in <strong>years BP</strong>. Select the concentration unit your data is in — the unit label is used for parameter hints only; data enters the model unconverted. Isotope data is <strong>optional</strong>.
      </div>

      <!-- Depth / Age -->
      <div class="card">
        <div class="card-title">📅 Depth / Age Model</div>
        <div class="upload-grid">
          <div>
            <div class="dropzone" id="dz-depth_age" ondragover="dzDrag(event,this)" ondragleave="dzLeave(this)" ondrop="dzDrop(event,'depth_age',this)">
              <input type="file" accept=".csv" onchange="dzFile(event,'depth_age',this.parentElement)">
              <div class="dz-icon">🗓</div>
              <div class="dz-label">Depth / Age</div>
              <div class="dz-hint">click or drag CSV</div>
              <div class="dz-name" id="name-depth_age"></div>
            </div>
            <div class="col-selector" id="cols-depth_age">
              <div style="display:grid;grid-template-columns:1fr 90px;gap:8px;align-items:end">
                <div>
                  <label>Depth column</label>
                  <select id="sel-depth_age-depth" onchange="storeCol('depth_age','depth',this.value);renderAgePlot()"></select>
                </div>
                <div>
                  <label>Unit</label>
                  <select id="da-depth-unit-sel" style="width:100%" onchange="renderAgePlot()">
                    <option value="cm" selected>cm</option>
                    <option value="mm">mm</option>
                  </select>
                </div>
              </div>
              <label>Age column (yrs BP)</label>
              <select id="sel-depth_age-age" onchange="storeCol('depth_age','age',this.value);renderAgePlot()"></select>
              <label>Age error column (2σ)</label>
              <select id="sel-depth_age-age_err" onchange="storeCol('depth_age','age_err',this.value)"></select>
            </div>
          </div>
        </div>
        <!-- Age-depth plot -->
        <div id="age-plot-wrap" style="display:none;margin-top:16px">
          <div style="font-size:11px;color:var(--muted);margin-bottom:6px;display:flex;justify-content:space-between;align-items:center;gap:10px">
            <span>Age–depth model &amp; extrapolation</span>
            <div style="display:flex;align-items:center;gap:6px;flex-shrink:0">
              <label style="font-size:10px;color:var(--muted)">Fit:</label>
              <select id="age-fit-type" onchange="renderAgePlot()"
                      style="font-size:10px;background:var(--surface);color:var(--text);
                             border:1px solid var(--border);border-radius:4px;padding:2px 6px">
                <option value="pchip" selected>Monotone spline (PCHIP)</option>
                <option value="linear">Linear regression</option>
              </select>
            </div>
            <span style="color:var(--accent);font-size:10px">Age range auto-populated ↓</span>
          </div>
          <div class="chart-wrap" style="height:260px;padding:10px">
            <canvas id="agePlotCanvas"></canvas>
          </div>
          <div style="margin-top:10px;display:grid;grid-template-columns:1fr 1fr;gap:10px">
            <div class="field">
              <label style="font-size:11px">Extrapolate min age (yrs BP)</label>
              <input type="number" id="age-extrap-min" style="width:100%"
                     oninput="onExtrapChange()" placeholder="auto">
            </div>
            <div class="field">
              <label style="font-size:11px">Extrapolate max age (yrs BP)</label>
              <input type="number" id="age-extrap-max" style="width:100%"
                     oninput="onExtrapChange()" placeholder="auto">
            </div>
          </div>
          <p style="font-size:10px;color:var(--muted);margin-top:6px;line-height:1.6">
            Fit line passes through all dated points. Dashed extensions show extrapolation beyond
            the data range using the end tangent. Edit the min/max fields or the extrap inputs to
            override the auto-populated calibration age range.
          </p>

          <!-- Growth rate chart + hiatus detection -->
          <div id="growth-rate-panel" style="margin-top:16px;padding-top:14px;border-top:1px solid var(--border)">
            <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:8px">
              <span style="font-size:12px;font-weight:600;color:var(--texthi)">📏 Growth Rate &amp; Hiatus Detection</span>
            </div>
            <div class="chart-wrap" style="height:160px;padding:10px">
              <canvas id="growthRateCanvas"></canvas>
            </div>
            <div style="display:grid;grid-template-columns:1fr 1fr;gap:10px;margin-top:8px;align-items:end">
              <div class="field">
                <label>Hiatus threshold (mm/yr)</label>
                <input type="number" id="hiatus-threshold" value="0.005" step="0.001" min="0"
                       style="width:100%" oninput="updateGrowthRate()">
                <div style="font-size:9.5px;color:var(--muted);margin-top:2px">
                  Intervals below this growth rate are flagged as potential hiatuses.
                </div>
              </div>
              <button onclick="autoDetectHiatuses()" class="btn btn-primary"
                      style="padding:8px 14px;font-size:11px;height:fit-content">Auto-detect hiatuses</button>
            </div>
            <div id="hiatus-zones" style="margin-top:12px">
              <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:6px">
                <span style="font-size:11px;font-weight:600;color:var(--texthi)">Exclusion zones</span>
                <button onclick="addHiatusZone()" class="btn btn-ghost" style="padding:3px 10px;font-size:11px">+ Add zone</button>
              </div>
              <div id="hiatus-zone-list"></div>
              <div id="hiatus-zone-hint" style="font-size:10px;color:var(--muted);margin-top:4px">
                No exclusion zones defined. Use auto-detect or add manually.
              </div>
            </div>
          </div>
        </div>
      </div>

      <!-- Trace Elements -->
      <div class="card">
        <div class="card-title">🔬 Trace Elements</div>
        <div class="dropzone" id="dz-trace_elem1" ondragover="dzDrag(event,this)" ondragleave="dzLeave(this)" ondrop="dzDrop(event,'trace_elem1',this)">
          <input type="file" accept=".csv" onchange="dzFile(event,'trace_elem1',this.parentElement)">
          <div class="dz-icon">⚗️</div>
          <div class="dz-label">Trace Elements CSV</div>
          <div class="dz-hint">click or drag — one CSV with depth + one or more proxy columns</div>
          <div class="dz-name" id="name-trace_elem1"></div>
        </div>
        <div id="te-config" style="display:none; margin-top:12px">
          <div style="display:grid;grid-template-columns:1fr 120px;gap:8px;margin-bottom:10px;align-items:end">
            <div>
              <label style="font-size:11px;color:var(--muted)">Depth column</label>
              <select id="te-depth-sel" style="width:100%;margin-top:4px" onchange="teDepthChanged(this.value)"></select>
            </div>
            <div>
              <label style="font-size:11px;color:var(--muted)">Depth unit</label>
              <select id="te-depth-unit-sel" style="width:100%;margin-top:4px" onchange="teDepthUnitChanged(this.value)">
                <option value="mm">mm</option>
                <option value="cm" selected>cm</option>
              </select>
            </div>
          </div>
          <div style="font-size:11px;color:var(--muted);margin-bottom:8px;padding-bottom:8px;border-bottom:1px solid var(--border)">
            Proxy columns — add one row per trace element:
          </div>
          <div id="te-rows"></div>
          <button onclick="addTERow()" class="btn btn-ghost" style="margin-top:4px;font-size:11px;padding:5px 14px">+ Add trace element</button>

          <!-- Preprocessing panel -->
          <div id="te-preproc" style="display:none;margin-top:16px;padding-top:14px;border-top:1px solid var(--border)">
            <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:10px">
              <span style="font-size:11px;font-weight:600;color:var(--texthi)">🔧 Data Preprocessing</span>
              <span id="te-row-badge" style="font-size:10px;color:var(--muted)"></span>
            </div>
            <div style="height:220px;margin-bottom:10px">
              <canvas id="te-preview-chart" style="width:100%;height:100%"></canvas>
            </div>
            <div style="display:grid;grid-template-columns:1fr 1fr 1fr 120px;gap:10px;align-items:end">
              <div class="field">
                <label style="font-size:11px">Block-average to N rows</label>
                <input type="number" id="te-target-n" min="10" step="10"
                       style="width:100%" oninput="updatePreviewChart()">
              </div>
              <div class="field">
                <label style="font-size:11px">Sigma-clip threshold (σ)</label>
                <input type="number" id="te-sigma" value="3" min="1" max="10" step="0.5"
                       style="width:100%" oninput="updatePreviewChart()">
              </div>
              <div class="field">
                <label style="font-size:11px">Window size (pts)
                  <span style="color:var(--muted);font-weight:400" title="Number of points in the centred rolling window used to compute local mean and SD for spike detection. Use 0 for global sigma-clip.">(?)</span>
                </label>
                <input type="number" id="te-winsize" value="50" min="0" step="10"
                       style="width:100%" oninput="updatePreviewChart()">
              </div>
              <button onclick="applyPreprocessing()" class="btn btn-primary"
                      style="padding:8px 14px;font-size:11px;height:fit-content">Apply &amp; Save</button>
            </div>
            <div id="te-preproc-status" style="font-size:10px;color:var(--muted);margin-top:6px;min-height:14px"></div>
          </div>

          <!-- TE cross-correlation scatter -->
          <div id="te-scatter-panel" style="display:none;margin-top:16px;padding-top:14px;border-top:1px solid var(--border)">
            <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:8px">
              <span style="font-size:12px;font-weight:600;color:var(--texthi)">📊 TE Cross-Correlation</span>
              <button onclick="document.getElementById('te-scatter-help').style.display=document.getElementById('te-scatter-help').style.display==='none'?'':'none'"
                      class="btn btn-ghost" style="padding:3px 8px;font-size:11px">?</button>
            </div>
            <div id="te-scatter-help" style="display:none;margin-bottom:10px;padding:8px 10px;
                 background:var(--bg);border:1px solid var(--border);border-radius:6px;font-size:11px;color:var(--text);line-height:1.7">
              If two trace elements are both controlled by drip rate via OMC dissociation
              kinetics, their speleothem concentrations should be correlated — higher drip
              rates reduce both. Positive correlation supports a shared drip rate control;
              low or no correlation suggests different processes dominate (redox, PCP, etc.).
              Use this to screen element pairs before multi-TE reconstruction.
            </div>
            <div style="display:grid;grid-template-columns:1fr 1fr 1fr;gap:8px;margin-bottom:8px;align-items:end">
              <div class="field">
                <label>X axis</label>
                <select id="te-scatter-x" onchange="updateScatterPlot()" style="width:100%"></select>
              </div>
              <div class="field">
                <label>Y axis</label>
                <select id="te-scatter-y" onchange="updateScatterPlot()" style="width:100%"></select>
              </div>
              <div class="field">
                <label>Fit</label>
                <select id="te-scatter-fit" onchange="updateScatterPlot()" style="width:100%">
                  <option value="none">None</option>
                  <option value="linear" selected>Linear</option>
                  <option value="exp">Exponential</option>
                </select>
              </div>
            </div>
            <div style="height:220px">
              <canvas id="te-scatter-chart"></canvas>
            </div>
            <div id="te-scatter-stats" style="font-size:10px;color:var(--muted);margin-top:4px"></div>
          </div>
        </div>
      </div>

      <!-- Isotopes -->
      <div class="card">
        <div class="card-title">💧 Isotope Data <span style="font-size:11px;font-weight:400;color:var(--muted)">(optional)</span></div>
        <div class="upload-grid">
          <div>
            <div class="dropzone" id="dz-isotope1" ondragover="dzDrag(event,this)" ondragleave="dzLeave(this)" ondrop="dzDrop(event,'isotope1',this)">
              <input type="file" accept=".csv" onchange="dzFile(event,'isotope1',this.parentElement)">
              <div class="dz-icon">💧</div>
              <div class="dz-label">Isotope (e.g. δ¹⁸O)</div>
              <div class="dz-hint">optional — click or drag CSV</div>
              <div class="dz-name" id="name-isotope1"></div>
            </div>
            <div class="col-selector" id="cols-isotope1">
              <label>Depth column</label>
              <select id="sel-isotope1-depth" onchange="storeCol('isotope1','depth',this.value)"></select>
              <label>Isotope column (‰)</label>
              <select id="sel-isotope1-proxy" onchange="storeCol('isotope1','proxy',this.value)"></select>
            </div>
          </div>
        </div>
      </div>

      <!-- Station settings -->
      <div class="card">
        <div class="card-title">🗺 Site Settings</div>
        <div class="form-grid">
          <div class="field">
            <label>Station name</label>
            <input type="text" id="station_name" value="Cave Site" placeholder="e.g. Heshang">
          </div>
          <div class="field">
            <label>Calibration age min (yrs BP)
              <span style="color:var(--muted);font-weight:400;font-size:10px"> — auto from plot</span>
            </label>
            <input type="number" id="calage_min" value="-50"
                   oninput="syncExtrapFromFields()">
          </div>
          <div class="field">
            <label>Calibration age max (yrs BP)
              <span style="color:var(--muted);font-weight:400;font-size:10px"> — auto from plot</span>
            </label>
            <input type="number" id="calage_max" value="9500"
                   oninput="syncExtrapFromFields()">
          </div>
        </div>
      </div>

      <div style="display:flex; gap:10px;">
        <button class="btn btn-primary" onclick="showPanel('params')">Next: Model Parameters →</button>
      </div>
    </div>


    <!-- MODEL PARAMETERS -->
    <div id="panel-params" class="panel">
      <div class="page-title">Model Parameters</div>
      <p class="page-desc">Geochemical and kinetic parameters for the drip rate reconstruction model. Defaults are from the reference dataset.</p>

      <!-- Analysis mode toggle -->
      <div class="card" style="margin-bottom:16px">
        <div class="card-title">⚙️ Analysis Mode</div>
        <div style="display:grid;grid-template-columns:1fr 1fr;gap:10px;margin-bottom:12px">
          <button id="mode-btn-full" onclick="setAnalysisMode('full')"
                  class="btn btn-primary"
                  style="text-align:left;padding:12px 14px;height:auto;display:flex;flex-direction:column;gap:4px">
            <span style="font-size:12px;font-weight:600">Full Quantification</span>
            <span style="font-size:10px;opacity:0.75;font-weight:400">Absolute drip rate (drips min⁻¹). Requires solution chemistry (aq. conc., Ca conc.).</span>
          </button>
          <button id="mode-btn-semi" onclick="setAnalysisMode('semi')"
                  class="btn btn-ghost"
                  style="text-align:left;padding:12px 14px;height:auto;display:flex;flex-direction:column;gap:4px">
            <span style="font-size:12px;font-weight:600">Semi-Quantitative</span>
            <span style="font-size:10px;opacity:0.75;font-weight:400">Relative drip rate (% of reference). No water chemistry needed.</span>
          </button>
        </div>
        <div id="mode-callout-full" style="font-size:11px;color:var(--muted);line-height:1.7">
          Solution concentrations and calibrated rate constants are required. Outputs absolute
          drip rate reconstructions suitable for hydrological quantification.
        </div>
        <div id="mode-callout-semi" style="display:none;font-size:11px;color:var(--muted);line-height:1.7">
          Element fractions and Kd values still control the <em>shape</em> of the response and should
          be set from literature. Solution chemistry inputs are hidden. Output is normalised to a
          reference value — either the record maximum or an optional user-supplied anchor drip rate
          (e.g. from modern monitoring).
          <div style="margin-top:10px;display:grid;grid-template-columns:1fr 1fr;gap:10px">
            <div class="field">
              <label>Anchor drip rate <span style="font-weight:400;color:var(--muted)">(optional, drips min⁻¹)</span></label>
              <input type="number" id="semi_anchor" placeholder="leave blank = normalise to max" min="0" step="0.1" style="width:100%">
            </div>
            <div class="field">
              <label>Reference period <span style="font-weight:400;color:var(--muted)">(optional, yrs BP)</span></label>
              <div style="display:flex;gap:6px">
                <input type="number" id="semi_ref_min" placeholder="from" style="flex:1">
                <input type="number" id="semi_ref_max" placeholder="to"   style="flex:1">
              </div>
            </div>
          </div>
        </div>
      </div>

      <!-- Element guidance note -->
      <div class="callout" style="margin-bottom:16px; font-size:11px; line-height:1.75">
        <strong style="color:var(--accent)">Element selection &amp; fraction guidance</strong><br>
        Elements exhibiting dissociation kinetics with organic metal complexes (OMC) are primarily
        d-block transition metals. However, redox-active elements (Mn, Fe) are not suitable as their
        speciation is governed by redox rather than drip rate. Elements beyond Ni, Co and Cu are not
        yet characterised in cave systems, and their partitioning behaviour is largely unexplored or
        inhibits utility as drip rate proxies — for example, Cu²⁺ has a high K<sub>d</sub> that
        limits sensitivity. Kp = −1 instructs the model to calculate Kp theoretically from cave temperature
        using the Wang &amp; Xu (2001) lattice strain model.<br><br>
        <strong style="color:var(--accent2)">Solution fractions:</strong>
        Total aqueous metal is partitioned into three components —
        <strong>X<sub>I</sub></strong> (inert: not bioavailable or complexed, does not participate in OMC dissociation),
        <strong>X<sub>F</sub></strong> (fast: labile but dissociates too rapidly to fractionate with drip rate), and
        <strong>X<sub>L</sub></strong> (labile: undergoes drip-rate-sensitive OMC dissociation; auto-calculated as
        1 − X<sub>I</sub> − X<sub>F</sub>). For Ni the recommended defaults are X<sub>I</sub> = 0.10,
        X<sub>F</sub> = 0.01, X<sub>L</sub> = 0.89.
      </div>

      <!-- Cave conditions -->
      <div class="card">
        <div class="card-title">🌡 Cave Conditions</div>
        <div class="form-grid">
          <div class="field">
            <label>Cave temperature (°C)</label>
            <input type="number" id="temp_C" value="16.5" step="0.1">
          </div>
          <div class="field">
            <label>Monitored drip rate (drips min⁻¹)
              <span style="color:var(--muted);font-weight:400;font-size:9px;cursor:help"
                    title="Modern observed drip rate used to anchor Kd estimation for all TEs. Leave blank or 0 to use literature Kd values."> (?)</span>
            </label>
            <input type="number" id="global_drip_rate" value="10" step="0.5" min="0"
                   oninput="updateAllParamHints()">
            <div id="global_drip_hint" style="font-size:9.5px;color:var(--muted);margin-top:3px;min-height:13px;line-height:1.4">
              Used to compute ln(Kd) from observed TE concentrations via kinetic inversion.
            </div>
          </div>
          <div class="field fullonly" style="grid-column:1/-1">
            <label style="font-size:11px;font-weight:600;color:var(--texthi);margin-bottom:6px">
              Ca aqueous concentration</label>
            <div style="display:flex;gap:6px;margin-bottom:8px">
              <button id="ca-mode-manual" class="mode-btn active"
                      onclick="setCaMode('manual')">Manual entry</button>
              <button id="ca-mode-csv" class="mode-btn"
                      onclick="setCaMode('csv')">Upload monitoring CSV</button>
            </div>
            <!-- Manual entry path (default) -->
            <div id="ca-manual-block">
              <div style="display:grid;grid-template-columns:1fr 1fr 1fr;gap:8px">
                <div>
                  <label style="font-size:10px;color:var(--muted)">Mean [Ca_aq]</label>
                  <input type="number" id="ca_conc" value="66.46" step="0.1"
                         oninput="updateCaHint();updateAllParamHints();fitCaPriorFromManual()">
                </div>
                <div>
                  <label style="font-size:10px;color:var(--muted)">Std dev (optional)</label>
                  <input type="number" id="ca_conc_sd" placeholder="blank = fixed"
                         step="0.1" oninput="fitCaPriorFromManual()">
                </div>
                <div>
                  <label style="font-size:10px;color:var(--muted)">Unit</label>
                  <select id="ca_unit"
                          onchange="updateCaHint();updateAllParamHints();fitCaPriorFromManual()">
                    <option value="ppb">ppb (µg/L)</option>
                    <option value="ppm" selected>ppm (mg/L)</option>
                  </select>
                </div>
              </div>
            </div>
            <!-- CSV upload path -->
            <div id="ca-csv-block" style="display:none">
              <div class="dropzone" id="dz-ca_aq"
                   ondragover="dzDrag(event,this)" ondragleave="dzLeave(this)"
                   ondrop="dzDrop(event,'ca_aq',this)">
                <input type="file" accept=".csv"
                       onchange="dzFile(event,'ca_aq',this.parentElement)">
                <div class="dz-label">Drop Ca monitoring CSV</div>
                <div class="dz-name" id="name-ca_aq"></div>
              </div>
              <div id="ca-col-selector" style="display:none;margin-top:8px">
                <label style="font-size:10px;color:var(--muted)">Ca column</label>
                <select id="ca_csv_col" onchange="fitCaPriorFromCsv()"
                        style="width:100%;margin-bottom:4px"></select>
                <label style="font-size:10px;color:var(--muted)">Unit in file</label>
                <select id="ca_csv_unit" onchange="fitCaPriorFromCsv()" style="width:100%">
                  <option value="ppb">ppb (µg/L)</option>
                  <option value="ppm" selected>ppm (mg/L)</option>
                </select>
              </div>
              <!-- Distribution chart -->
              <div style="display:flex;align-items:center;justify-content:space-between;margin-top:10px;margin-bottom:4px">
                <span style="font-size:10px;color:var(--muted)">Distribution fit</span>
                <button id="ca-dist-chart_toggle" onclick="toggleLogScale('ca-dist-chart')"
                        style="background:none;border:1px solid var(--border);color:var(--muted);
                               border-radius:4px;cursor:pointer;font-size:10px;padding:2px 8px">Log x</button>
              </div>
              <div style="height:160px;display:none" id="ca-dist-chart-wrap">
                <canvas id="ca-dist-chart"></canvas>
              </div>
              <div id="ca-dist-desc" style="display:none;margin-top:8px;padding:8px 10px;
                   background:var(--bg);border:1px solid var(--border);border-radius:6px;font-size:10px"></div>
            </div>
            <!-- Prior summary -->
            <div id="ca-prior-summary" style="font-size:10px;color:var(--muted);margin-top:8px;
                 padding:6px 8px;background:var(--bg);border:1px solid var(--border);
                 border-radius:4px;display:none"></div>
            <div id="ca_hint" style="font-size:9.5px;color:var(--muted);margin-top:3px;min-height:13px;line-height:1.4"></div>
          </div>
        </div>
        <div class="callout" style="margin-top:12px;font-size:11px;line-height:1.75">
          <strong style="color:var(--accent)">K<sub>d</sub> estimation</strong><br>
          Each trace element card below offers two modes for setting ln(K<sub>d</sub>):<br>
          <strong>From drip rate</strong> — for monitored sites: enter the measured drip rate and K<sub>d</sub>
          is back-calculated from the observed calcite concentration using OMC dissociation kinetics
          (<code style="display:inline;padding:2px 4px;margin:0">f = 1 − e<sup>−K<sub>e</sub>·τ</sup></code>,
          τ = 60/V).<br>
          <strong>Literature / manual</strong> — for unmonitored sites: enter ln(K<sub>d</sub>) directly
          from literature values or your own calibration.
        </div>

      <!-- TE parameter cards (dynamically generated) -->
      <div id="te-all-param-cards"></div>

      </div>

      <div style="display:flex; gap:10px;">
        <button class="btn btn-ghost" onclick="showPanel('data')">← Back</button>
        <button class="btn btn-primary" onclick="showPanel('output')">Next: Output Options →</button>
      </div>
    </div>


    <!-- OUTPUT OPTIONS -->
    <div id="panel-output" class="panel">
      <div class="page-title">Output Options</div>
      <p class="page-desc">Configure what gets saved and how the realisations are generated.</p>

      <div class="card">
        <div class="card-title">📊 Realisations (for RQA)</div>
        <div class="field" style="margin-bottom:12px">
          <label style="display:flex;align-items:center;gap:10px;cursor:pointer;color:var(--text)">
            <input type="checkbox" id="generate_realisations"
                   style="accent-color:var(--accent);width:14px;height:14px"
                   onchange="document.getElementById('rqa-options').style.display=this.checked?'':' none'">
            Generate realisations CSV (needed for RQA)
          </label>
        </div>
        <div id="rqa-options">
        <div class="form-grid">
          <div class="field">
            <label>Number of realisations</label>
            <input type="number" id="n_realisations" value="1000" min="10" max="10000" step="10">
          </div>
          <div class="field">
            <label>Random seed (reproducibility)</label>
            <input type="number" id="rng_seed" value="42">
          </div>
        </div>
        </div>
        <div style="margin-top:12px; font-size:11px; color:var(--muted); line-height:1.7">
          Each realisation is an independent draw from the joint posterior PDF of drip rate at each time step,
          forming a complete plausible time series. These are saved as a CSV with shape
          <em>(N_timesteps × N_realisations + 1)</em> suitable for direct input to RQA tools.
        </div>
      </div>

      <div class="card">
        <div class="card-title">📋 Summary Statistics</div>
        <div style="font-size:12px; color:var(--text); line-height:1.8">
          The summary CSV is always saved and contains the following percentile columns for each time step:
          <br><br>
          <code style="color:var(--accent)">age, pc05, pc10, pc25, pc50, pc75, pc90, pc95</code>
        </div>
      </div>

      <div class="card">
        <div class="card-title">⚡ Proxy Record</div>
        <div class="field" style="margin-bottom:12px">
          <label style="display:flex; align-items:center; gap:10px; cursor:pointer; color:var(--text)">
            <input type="checkbox" id="use_cached_proxy" style="accent-color:var(--accent); width:14px; height:14px;">
            Use cached proxy record (ProxyRecord.pkl) if available
          </label>
        </div>
        <div class="callout warn" style="margin:0">
          The proxy record (BayProX) is the slow step. Progress and ETA are shown on the Run page.
          If you have already run it once and are only changing model parameters (Kd, Kp, etc.),
          check this box to skip straight to the drip rate step.
        </div>
      </div>

      <div class="card">
        <div class="card-title">🔧 Advanced — Drip Rate Grid</div>
        <div style="font-size:11px;color:var(--muted);line-height:1.7;margin-bottom:12px">
          The forward model evaluates h(V) on a grid from V<sub>min</sub> to V<sub>max</sub>.
          If the model's predicted [TE]<sub>calcite</sub> at V<sub>max</sub> is still above
          the observed data range, increase V<sub>max</sub>. This is needed when K₀ is large
          (high Kp or high aq/ca ratio) and most of the kinetic fractionation occurs at
          high drip rates beyond the default 100 drips/min.
        </div>
        <div class="form-grid">
          <div class="field">
            <label>V<sub>max</sub> (drips min⁻¹)</label>
            <div style="display:flex;gap:6px;align-items:center">
              <input type="number" id="v_max" value="100" min="10" max="10000" step="10" style="flex:1">
              <button class="btn btn-ghost" onclick="autoSetVmax()" style="padding:5px 10px;font-size:11px"
                      title="Set VMAX based on where h(V) crosses the data range">Auto</button>
            </div>
            <div id="v_max_hint" style="font-size:9.5px;color:var(--muted);margin-top:3px;min-height:13px;line-height:1.4">
              Default: 100. Use Auto to fit to data range.
            </div>
          </div>
          <div class="field">
            <label>V grid resolution</label>
            <input type="number" id="v_res" value="5000" min="500" max="50000" step="500">
            <div style="font-size:9.5px;color:var(--muted);margin-top:3px">
              Default: 5000. Higher = finer grid but slower.
            </div>
          </div>
        </div>
      </div>

      <div style="display:flex; gap:10px;">
        <button class="btn btn-ghost" onclick="showPanel('params')">← Back</button>
        <button class="btn btn-primary" onclick="showPanel('run')">Next: Run Model →</button>
      </div>
    </div>


    <!-- RUN -->
    <div id="panel-run" class="panel">
      <div class="page-title">Run Model</div>
      <p class="page-desc">Review your configuration and launch the reconstruction pipeline.</p>

      <div class="card" id="configSummary">
        <div class="card-title">📋 Configuration Summary</div>
        <div id="summaryContent" style="font-size:12px; line-height:1.9; color:var(--text)"></div>
      </div>

      <div class="card" id="runtime-estimate-card">
        <div class="card-title">⏱ Runtime Estimate</div>
        <div id="runtime-estimate-body" style="font-size:12px;color:var(--text);line-height:1.9"></div>
      </div>

      <div class="card" id="sanity-card" style="display:none">
        <div class="card-title">⚠️ Parameter Warnings</div>
        <div id="sanity-body" style="font-size:11px;line-height:1.8"></div>
      </div>

      <div class="card">
        <div class="card-title">Progress</div>
        <div id="stageLabel" style="font-size:12px; color:var(--text); min-height:18px"></div>
        <div class="progress-wrap">
          <div class="progress-bar" id="progressBar"></div>
        </div>
        <div class="progress-label">
          <span id="progressPct">0%</span>
          <span id="etaLabel"></span>
        </div>
        <hr class="divider">
        <div class="card-title">Log</div>
        <div class="log-box" id="logBox"></div>
      </div>

      <div style="display:flex; gap:10px; flex-wrap:wrap;">
        <button class="btn btn-primary" id="runBtn" onclick="startRun()">▶ Start Run</button>
        <button class="btn btn-ghost" onclick="showPanel('output')">← Back</button>
        <button class="btn btn-ghost" id="resultsBtn" style="display:none" onclick="showPanel('results')">
          View Results →
        </button>
      </div>
    </div>


    <!-- RESULTS -->
    <div id="panel-results" class="panel">
      <div class="page-title">Output Explorer</div>
      <p class="page-desc">Reconstructed drip rate and hydrological classification.</p>

      <!-- Tab bar -->
      <div style="display:flex;gap:0;border-bottom:1px solid var(--border);margin-bottom:16px">
        <button id="tab-ts" onclick="switchResultTab('ts')" class="result-tab active"
                style="padding:8px 18px;font-size:11px;background:none;border:none;
                       border-bottom:2px solid var(--accent);color:var(--texthi);cursor:pointer">📈 Time Series</button>
        <button id="tab-sf" onclick="switchResultTab('sf')" class="result-tab"
                style="padding:8px 18px;font-size:11px;background:none;border:none;
                       border-bottom:2px solid transparent;color:var(--muted);cursor:pointer">💧 Smart &amp; Friedrich</button>
        <button id="tab-pdf" onclick="switchResultTab('pdf')" class="result-tab"
                style="padding:8px 18px;font-size:11px;background:none;border:none;
                       border-bottom:2px solid transparent;color:var(--muted);cursor:pointer">🌈 PDF Heatmap</button>
        <button id="tab-age" onclick="switchResultTab('age')" class="result-tab"
                style="padding:8px 18px;font-size:11px;background:none;border:none;
                       border-bottom:2px solid transparent;color:var(--muted);cursor:pointer">📐 Age Model</button>
      </div>

      <!-- Y-axis mode toggle -->
      <div style="display:flex;align-items:center;gap:10px;margin-bottom:14px;padding:6px 10px;
                  background:var(--surface);border:1px solid var(--border);border-radius:6px">
        <label style="font-size:11px;color:var(--texthi);font-weight:600">Y-axis:</label>
        <button id="ymode-drip" class="mode-btn active" onclick="setYMode('drip')">Drip rate (min⁻¹)</button>
        <button id="ymode-tau" class="mode-btn" onclick="setYMode('tau')">τ (s)</button>
        <span style="color:var(--muted);font-size:9px;cursor:help;margin-left:4px"
              title="Drip rates assume replacement of the thin-film by each subsequent drop on a stalagmite apex. For stalactites, flowstones or sub-aqueous samples, use τ (s), where τ = 60/d. τ represents the thin-film residence time and generalises the kinetic window to deposits where drip-replenishment doesn't cleanly define the accumulation interval.">(?)</span>
      </div>

      <!-- Time series panel -->
      <div id="rtab-ts">
        <div class="card">
          <div class="card-title" id="ts-card-title">📈 Drip Rate vs Time</div>
          <div class="chart-wrap">
            <canvas id="dripChart"></canvas>
          </div>
        </div>
      </div>

      <!-- Smart & Friedrich panel -->
      <div id="rtab-sf" style="display:none">
        <div class="card">
          <div class="card-title">💧 Smart &amp; Friedrich (1987) Classification</div>
          <div style="font-size:11px;color:var(--muted);margin-bottom:10px;line-height:1.7">
            Hydrological classification of cave drip sites by mean discharge and
            coefficient of variation (CV = σ/μ). Point shows median reconstruction;
            bars span the 25th–75th percentile range of the time series.
          </div>
          <div class="chart-wrap" style="height:380px">
            <canvas id="sfChart"></canvas>
          </div>
          <div id="sf-label" style="text-align:center;margin-top:10px;font-size:12px;
                                    color:var(--accent);font-weight:500"></div>
        </div>
        <div class="card" style="margin-top:12px">
          <div class="card-title" style="font-size:11px">Zone definitions</div>
          <div style="display:grid;grid-template-columns:1fr 1fr;gap:8px;font-size:11px;color:var(--text);line-height:1.7">
            <div><span style="color:#5b9bd5">■</span> <strong>Seepage / percolation</strong><br>
              Low mean discharge, low CV. Slow matrix flow through fine pores.
              Highly attenuated signal.</div>
            <div><span style="color:#70ad47">■</span> <strong>Fracture / conduit</strong><br>
              Low–moderate discharge, high CV. Episodic fracture recharge.
              Strong seasonal signal.</div>
            <div style="margin-top:6px"><span style="color:#ed7d31">■</span> <strong>Buffered overflow</strong><br>
              High discharge, low CV. Perched water table overflow.
              Sustained but seasonally modulated.</div>
            <div style="margin-top:6px"><span style="color:#a855f7">■</span> <strong>Flood / conduit overflow</strong><br>
              High discharge, high CV. Rapid conduit response to recharge events.</div>
          </div>
        </div>
      </div>

      <!-- PDF heatmap panel -->
      <div id="rtab-pdf" style="display:none">
        <div class="card">
          <div class="card-title">🌈 Drip Rate Probability Density</div>
          <div style="display:flex;gap:10px;margin-bottom:10px;align-items:center;flex-wrap:wrap">
            <label style="font-size:10px;color:var(--muted)">Colour map:</label>
            <select id="pdf-cmap" onchange="renderPdfHeatmap()"
                    style="font-size:11px;padding:3px 8px">
              <optgroup label="Sequential">
                <option value="greens">Greens</option>
                <option value="blues">Blues</option>
                <option value="reds">Reds</option>
              </optgroup>
              <optgroup label="Perceptual">
                <option value="viridis" selected>Viridis</option>
                <option value="inferno">Inferno</option>
                <option value="plasma">Plasma</option>
                <option value="magma">Magma</option>
                <option value="cividis">Cividis</option>
                <option value="turbo">Turbo</option>
              </optgroup>
            </select>
            <label style="font-size:10px;color:var(--muted);margin-left:10px">Log density:</label>
            <input type="checkbox" id="pdf-log" onchange="renderPdfHeatmap()"
                   style="accent-color:var(--accent)">
            <label style="font-size:10px;color:var(--muted);margin-left:10px">Percentiles:</label>
            <input type="checkbox" id="pdf-overlay" onchange="renderPdfHeatmap()"
                   style="accent-color:var(--accent)" checked>
          </div>
          <div style="height:400px;position:relative">
            <canvas id="pdfHeatmapCanvas" style="width:100%;height:100%"></canvas>
          </div>
          <div style="display:flex;justify-content:space-between;font-size:10px;color:var(--muted);margin-top:4px">
            <span id="pdf-xrange"></span>
            <span id="pdf-yrange"></span>
          </div>
        </div>
      </div>

      <!-- Age model panel -->
      <div id="rtab-age" style="display:none">
        <div class="card">
          <div class="card-title">📐 BayProX Age Model</div>
          <div style="font-size:11px;color:var(--muted);margin-bottom:10px;line-height:1.7">
            The age-depth relationship computed by BayProX from the U-Th dated
            points. Uncertainty envelope shows the model posterior spread.
          </div>
          <div class="chart-wrap">
            <canvas id="ageModelChart"></canvas>
          </div>
        </div>
      </div>

      <div style="display:flex;gap:10px;margin-top:4px">
        <button class="btn btn-ghost" onclick="showPanel('downloads')">⬇ Downloads</button>
      </div>
    </div>


    <!-- DOWNLOADS -->
    <div id="panel-downloads" class="panel">
      <div class="page-title">Downloads</div>
      <p class="page-desc">Output files from the most recent run.</p>

      <div class="output-list" id="downloadList">
        <div style="color:var(--muted); font-size:12px;">No outputs yet — run the model first.</div>
      </div>

    </div><!-- /panel-downloads -->


    <!-- ABOUT -->
    <div id="panel-about" class="panel">
      <div class="page-title">About</div>
      <p class="page-desc">Project background, authorship and funding acknowledgements.</p>

      <div class="card" style="border-color:var(--accent);border-width:2px">
        <div class="card-title" style="font-size:13px">🌍 Parent Project — QUEST</div>
        <div style="font-size:12px; line-height:1.9; color:var(--text)">
          <div style="display:flex;gap:16px;align-items:center;margin-bottom:14px;flex-wrap:wrap">
            <a href="https://quest.pik-potsdam.de/" target="_blank">
              <img src="https://quest.pik-potsdam.de/wp-content/uploads/2017/03/Logo-AI_QUEST-1024x527-300x154.jpg"
                   alt="QUEST — QUantitative palaeoEnvironments from SpeleoThems"
                   style="height:60px;border-radius:4px;background:#fff;padding:4px"
                   onerror="this.style.display='none'">
            </a>
            <div>
              <div style="font-family:var(--serif);font-size:18px;font-weight:600;color:var(--texthi)">
                QUEST</div>
              <div style="color:var(--accent);font-size:12px">
                QUantitative palaeoEnvironments from SpeleoThems</div>
            </div>
          </div>
          <div style="margin-bottom:12px">
            This tool was developed within the QUEST project — an international research
            collaboration that develops new techniques for extracting quantitative
            palaeoclimate information from speleothems. QUEST combines field and laboratory
            experiments on water/mineral chemistry with innovative physical and numerical
            analyses to deliver quantitative reconstructions of past hydrology and temperature.
          </div>
          <div style="display:flex;gap:12px;flex-wrap:wrap;align-items:center">
            <a href="https://quest.pik-potsdam.de/" target="_blank"
               class="btn btn-primary" style="font-size:11px;padding:6px 14px;text-decoration:none">
              Visit QUEST project ↗
            </a>
            <a href="https://cordis.europa.eu/project/id/691037" target="_blank"
               style="color:var(--accent);font-size:11px;text-decoration:none">
              EU CORDIS project record ↗
            </a>
          </div>
        </div>
      </div>

      <div class="card">
        <div class="card-title">📄 Manuscript</div>
        <div style="font-size:12px; line-height:1.9; color:var(--text)">
          <div style="font-size:14px; color:var(--texthi); margin-bottom:10px; line-height:1.5">
            Quantitative Holocene precipitation reconstruction from stalagmite trace metal kinetics reveals East Asian monsoon drivers
          </div>
          <div style="color:var(--muted); margin-bottom:12px">
            Adam Hartland, Bedartha Goswami, Jungho Park, Sebastian Höpker, Dorisel Torres Rojas,
            Jin Liao, Bethany R.S. Fox, Norbert Marwan, Sebastian F.M. Breitenbach, Chaoyong Hu
          </div>
          <div style="margin-bottom:8px">
            <em>Submitted to Nature Geoscience</em>
          </div>
          <div class="callout warn" style="margin-top:12px; font-size:11px">
            This is a pre-publication version of the drip rate model, provided for review purposes only.
            It should not be cited, distributed, or used for purposes other than peer review without
            permission from the corresponding author.
          </div>
          <div style="margin-top:14px; font-size:11px; color:var(--muted)">
            Correspondence: <a href="mailto:adam.hartland@lincolnagritech.co.nz"
              style="color:var(--accent)">adam.hartland@lincolnagritech.co.nz</a>
          </div>
          <div style="margin-top:8px; font-size:11px; color:var(--muted)">
            Code &amp; data:
            <a href="https://github.com/waikatosci/paleodriprates" target="_blank"
               style="color:var(--accent)">github.com/waikatosci/paleodriprates</a>
            &nbsp;·&nbsp;
            <a href="https://doi.org/10.5281/zenodo.16392750" target="_blank"
               style="color:var(--accent)">doi:10.5281/zenodo.16392750</a>
          </div>
        </div>
      </div>

      <div class="card">
        <div class="card-title">🔬 Project Description</div>
        <div style="font-size:12px; line-height:1.9; color:var(--text)">
          This tool implements a novel kinetic proxy approach that harnesses the dissociation
          of organic-metal complexes (OMC) in cave dripwater to reconstruct past drip rates —
          and thus precipitation — from stalagmite trace metals. The underlying model
          links the concentration of trace elements in speleothem calcite to past fluid residence
          times via the kinetics of metal-organic dissociation, calibrated through monitoring data
          and refined using Monte Carlo propagation of paleo-temperature uncertainty. While
          developed for speleothems, the approach generalises to any mineral deposit formed from
          waters containing OMC with sufficiently slow dissociation kinetics — see below.<br><br>
          The method was developed and validated using stalagmite HS4 from Heshang Cave, China,
          reconstructing Holocene rainfall for the Yangtze region and providing new insights
          into East Asian Summer Monsoon dynamics over the past 9,500 years. The probabilistic
          reconstruction framework — including the full ensemble of drip rate realisations
          available for download — enables recurrence quantification analysis (RQA) and other
          nonlinear time series methods.<br><br>
          <strong>Generalisation beyond stalagmite drip rates</strong><br>
          While the default output is expressed as drip rate (drips min⁻¹), the fundamental
          quantity recovered by the model is the thin-film residence time, τ = 60/<em>d</em> (s).
          For stalagmites, τ maps directly to drip rate because each drop replaces the thin
          film on the apex — but this is merely the simplest case of a general principle:
          wherever organic–metal complexes (OMC) are present in an aqueous system, their
          dissociation kinetics encode the timescale of fluid turnover.<br><br>

          <strong style="color:var(--accent)">Carbonate depositional settings</strong>
          <div style="margin:8px 0 8px 12px;padding:8px 12px;border-left:2px solid var(--border);font-size:11.5px;line-height:1.8">
            <strong>Stalagmites</strong> — τ = drip interval. Each drop replaces the thin film
            at the apex; the kinetic window is defined by drip-replenishment.<br>
            <strong>Stalactites</strong> — τ represents the fluid residence time as the film
            migrates along the soda-straw or curtain surface. The kinetic window is set by
            film thickness, flow path length, and gravity-driven velocity.<br>
            <strong>Flowstones</strong> — τ is the sheet-flow residence time across the mineral
            surface, governed by slope, roughness, and volumetric flux.<br>
            <strong>Sub-aqueous carbonates</strong> — τ represents the equilibration timescale
            of the thin boundary layer with bulk solution via diffusion, applicable to
            lacustrine tufas, pool deposits, and submarine cements.<br>
            <strong>Long integration periods</strong> — where the kinetic accumulation window
            <em>M</em> exceeds the reservoir mixing time, τ approximates bulk fluid residence
            time. Calibration to volumetric flow then requires reservoir volume:
            <em>Q</em> = <em>V</em><sub>reservoir</sub> / τ (m³ s⁻¹).
          </div>

          <strong style="color:var(--accent)">Beyond carbonates</strong>
          <div style="margin:8px 0 8px 12px;padding:8px 12px;border-left:2px solid var(--border);font-size:11.5px;line-height:1.8">
            Organic ligands capable of forming kinetically inert metal complexes are ubiquitous
            in Earth's aqueous compartments — rivers, lakes, estuaries, groundwater, soils,
            ocean margins, and hydrothermal systems. Any mineral phase that incorporates
            trace metals from such waters is, in principle, a recorder of fluid residence time
            via the same dissociation-kinetic mechanism. The only fundamental requirement is
            that the OMC dissociation timescale is sufficiently slow relative to the fluid
            turnover rate to produce a measurable kinetic signal — i.e., that the system is
            not at secular equilibrium with respect to metal release.<br><br>
            Even at secular equilibrium, the method resolves a <em>minimum</em> fluid residence
            time (or equivalently, a maximum flow rate): the observation that TE concentrations
            have reached the equilibrium ceiling K₀ constrains τ to be at least several
            multiples of the dissociation half-life. This provides a one-sided bound on
            paleohydrology that is still quantitatively useful — analogous to how radiocarbon
            provides a minimum age when a sample is beyond the dating range.
          </div>

          <strong style="color:var(--accent)">Toward absolute paleohydrology</strong>
          <div style="margin:8px 0 8px 12px;padding:8px 12px;border-left:2px solid var(--border);font-size:11.5px;line-height:1.8">
            Conventional paleohydrology relies on empirical transfer functions or semi-quantitative
            indices — proxies that track <em>relative</em> wetness or aridity without a
            mechanistic connection to physical flow rates. The kinetic approach here grounds
            hydrological reconstruction in the same class of absolute physical measurement
            that underpins geochronology: a rate law (OMC dissociation) with known kinetic
            parameters, operating over a measurable timescale (τ), recorded in a datable
            mineral archive. Just as radiometric decay provides an absolute clock, OMC
            dissociation kinetics provide an absolute flowmeter — extending paleohydrology
            from relative proxy interpretation into the methodological framework of
            geochronometry.
          </div>
          Use the τ (s) toggle in the Output Explorer to display results as residence time
          instead of drip rate.
        </div>
      </div>

      <div class="card">
        <div class="card-title">🏛 Affiliations</div>
        <div style="font-size:12px; line-height:1.9; color:var(--text)">
          <div style="display:grid; grid-template-columns: 1fr 1fr; gap:8px 24px;">
            <div>Lincoln Agritech Ltd, Ruakura, Hamilton, New Zealand</div>
            <div>Te Aka Mātuatua, School of Science, University of Waikato, New Zealand</div>
            <div>Lincoln University, Lincoln, New Zealand</div>
            <div>Potsdam Institute for Climate Impact Research, Potsdam, Germany</div>
            <div>IISER Pune, India</div>
            <div>China University of Geosciences, Wuhan, China</div>
            <div>University of Huddersfield, UK</div>
            <div>Northumbria University, Newcastle upon Tyne, UK</div>
          </div>
        </div>
      </div>

      <div class="card">
        <div class="card-title">💰 Funding</div>
        <div style="font-size:12px; line-height:1.9; color:var(--text)">
          This research was supported by the following funders:
          <div style="margin-top:12px; display:flex; flex-direction:column; gap:10px;">

            <!-- EU / QUEST -->
            <div style="display:flex; gap:14px; align-items:center; padding:12px 14px;
                        background:var(--bg); border:1px solid var(--border); border-radius:8px">
              <div style="display:flex;gap:10px;align-items:center;flex-shrink:0">
                <img src="https://quest.pik-potsdam.de/wp-content/uploads/2017/03/eu_flag-1-300x200.png"
                     alt="EU Flag" style="height:36px;border-radius:2px"
                     onerror="this.outerHTML='<span style=font-size:28px>🇪🇺</span>'">
                <img src="https://quest.pik-potsdam.de/wp-content/uploads/2017/03/Logo-AI_QUEST-1024x527-300x154.jpg"
                     alt="QUEST" style="height:32px;border-radius:3px;background:#fff;padding:2px"
                     onerror="this.style.display='none'">
              </div>
              <div>
                <div style="color:var(--texthi); font-weight:500">EU Horizon 2020 — Marie Skłodowska-Curie Actions</div>
                <div style="color:var(--muted); font-size:11px">Grant no. 691037 —
                  <a href="https://quest.pik-potsdam.de/" target="_blank"
                     style="color:var(--accent)">QUEST: QUantitative paleoEnvironments from SpeleoThems</a>
                </div>
              </div>
            </div>

            <!-- Te Apārangi -->
            <div style="display:flex; gap:14px; align-items:center; padding:12px 14px;
                        background:var(--bg); border:1px solid var(--border); border-radius:8px">
              <div style="flex-shrink:0">
                <img src="https://www.royalsociety.org.nz/assets/Uploads/Logo-Te-Aparangi-lockup__FillMaxWzIwMCwyMDBd.png"
                     alt="Royal Society Te Apārangi" style="height:40px"
                     onerror="this.outerHTML='<span style=font-size:28px>🇳🇿</span>'">
              </div>
              <div>
                <div style="color:var(--texthi); font-weight:500">
                  <a href="https://www.royalsociety.org.nz/" target="_blank"
                     style="color:var(--texthi);text-decoration:none">Te Apārangi — Royal Society of New Zealand</a>
                </div>
                <div style="color:var(--muted); font-size:11px">Grant RIS-UOW1501 &nbsp;·&nbsp; Rutherford Discovery Fellowship RDF-UOW1601</div>
              </div>
            </div>

            <!-- MBIE -->
            <div style="display:flex; gap:14px; align-items:center; padding:12px 14px;
                        background:var(--bg); border:1px solid var(--border); border-radius:8px">
              <div style="flex-shrink:0">
                <img src="https://www.mbie.govt.nz/dmsdocument/5765-ministry-of-business-innovation-and-employment-logo"
                     alt="MBIE" style="height:36px"
                     onerror="this.outerHTML='<span style=font-size:28px>🇳🇿</span>'">
              </div>
              <div>
                <div style="color:var(--texthi); font-weight:500">
                  <a href="https://www.mbie.govt.nz/" target="_blank"
                     style="color:var(--texthi);text-decoration:none">Ministry for Business, Innovation and Employment (MBIE)</a>
                </div>
                <div style="color:var(--muted); font-size:11px">Grant UOWX2102</div>
              </div>
            </div>

          </div>
        </div>
      </div>
    </div>

    </div>

  </main>

  <!-- Footer -->
  <div class="app-footer">
    <span class="footer-badge">⚠ Pre-publication review</span>
    <div class="footer-text">
      This is a pre-publication version of the drip rate model provided for review purposes and is based on the following manuscript submitted to <em>Nature Geoscience</em>:
      <strong> Quantitative Holocene precipitation reconstruction from stalagmite trace metal kinetics reveals East Asian monsoon drivers</strong>
      — Hartland et al.
      &nbsp;·&nbsp;
      <a href="#" onclick="showPanel('about'); return false;">About &amp; Funding ↗</a>
    </div>
    <div style="text-align:right;font-size:9px;color:var(--muted);margin-top:4px" id="build-stamp">Build v44</div>
  </div>

</div>

<script>
// ── State ────────────────────────────────────────────────────────────────────
const colMap = {};   // colMap[key][role] = selected column name
let pollTimer = null;
let chart = null;
let runStartTime = null;

// ── Navigation ───────────────────────────────────────────────────────────────
function showPanel(name) {
  document.querySelectorAll('.panel').forEach(p => p.classList.remove('active'));
  document.querySelectorAll('.nav-item').forEach(n => n.classList.remove('active'));
  document.getElementById('panel-' + name).classList.add('active');
  const _main = document.querySelector('main');
  if (_main) _main.scrollTo({top: 0, behavior: 'smooth'});
  // Find matching nav item safely without template-literal querySelector
  document.querySelectorAll('.nav-item').forEach(n => {
    const oc = n.getAttribute('onclick') || '';
    if (oc.includes("'" + name + "'")) n.classList.add('active');
  });
  if (name === 'run') { buildSummary(); updateRuntimeEstimate(); runSanityChecks(); }
  if (name === 'params') updateAllParamHints();
  if (name === 'results') loadChart();
  if (name === 'downloads') refreshDownloads();
}

// ── Upload / dropzone helpers ────────────────────────────────────────────────
function dzDrag(e, el) { e.preventDefault(); el.classList.add('drag'); }
function dzLeave(el)   { el.classList.remove('drag'); }

function dzDrop(e, key, el) {
  e.preventDefault();
  el.classList.remove('drag');
  const file = e.dataTransfer.files[0];
  if (file) uploadFile(key, file, el);
}

function dzFile(e, key, el) {
  const file = e.target.files[0];
  if (file) uploadFile(key, file, el);
}

function uploadFile(key, file, el) {
  const fd = new FormData();
  fd.append(key, file);
  fetch('/upload', { method: 'POST', body: fd })
    .then(r => r.json())
    .then(data => {
      if (data[key] && data[key].columns) {
        el.classList.add('ready');
        document.getElementById('name-' + key).textContent = file.name;
        if (key === 'depth_age' && data[key].data) {
          agePlotData = data[key].data;
        }
        if (key === 'trace_elem1' && data[key].data) {
          initPreprocessing(data[key].data,
            Object.keys(data[key].data)[0],  // placeholder — updated after depth sel
            data[key].rows);
        }
        // Ca monitoring CSV
        if (key === 'ca_aq' && data[key].data) {
          window.caAqCsvData = data[key].data;
          const sel = document.getElementById('ca_csv_col');
          if (sel) {
            sel.innerHTML = data[key].columns.map(c =>
              `<option value="${c}">${c}</option>`).join('');
            document.getElementById('ca-col-selector').style.display = '';
            fitCaPriorFromCsv();
          }
        }
        // TE monitoring CSVs (match te1_aq, te2_aq, etc.)
        const aqMatch = key.match(/^(te\d+)_aq$/);
        if (aqMatch && data[key].data) {
          const p = aqMatch[1];
          window[`${p}AqCsvData`] = data[key].data;
          const sel = document.getElementById(`${p}_aq_csv_col`);
          if (sel) {
            sel.innerHTML = data[key].columns.map(c =>
              `<option value="${c}">${c}</option>`).join('');
            document.getElementById(`${p}-aq-col-selector`).style.display = '';
            fitAqPriorFromCsv(p);
          }
        }
        populateColSelectors(key, data[key].columns);
        if (key === 'depth_age') {
          // Auto-detect depth unit from column name
          const daDepCol = document.getElementById('sel-depth_age-depth')?.value || '';
          const daSel = document.getElementById('da-depth-unit-sel');
          if (daSel) {
            const dcl = daDepCol.toLowerCase();
            if (dcl.includes('(mm)') || dcl.includes('_mm') || dcl.endsWith('mm') || dcl.includes(' mm')) {
              daSel.value = 'mm';
            } else if (dcl.includes('(cm)') || dcl.includes('_cm') || dcl.endsWith('cm') || dcl.includes(' cm')) {
              daSel.value = 'cm';
            }
          }
          renderAgePlot();
        }
      }
    })
    .catch(err => console.error('Upload error for', key, err));
}

// ── TE row state ─────────────────────────────────────────────────────────
let teColumns = [];   // all columns from uploaded TE CSV
let teDepthCol = '';  // selected depth column
let teDepthUnit = 'cm'; // auto-detected: 'mm' or 'cm'
let teRowData  = [{col: '', unit: 'ppm'}];  // [{col, unit}, ...] — default TE1

const UNIT_OPTS = [
  {v:'ppb', l:'ppb  (ng/g)'},
  {v:'ppm', l:'ppm  (µg/g)'},
  {v:'ug/g', l:'µg/g'},
  {v:'mg/kg', l:'mg/kg'},
];
const ELEM_OPTS = `
  <optgroup label="d-block (well characterised)">
    <option value="Ni">Ni — Nickel</option>
    <option value="Co">Co — Cobalt</option>
    <option value="Cu">Cu — Copper (high Kd)</option>
  </optgroup>
  <optgroup label="d-block (literature Kp)">
    <option value="Zn">Zn — Zinc</option>
    <option value="Mn">Mn — Manganese (redox-sensitive)</option>
    <option value="Fe">Fe — Iron (redox-sensitive)</option>
    <option value="Cd">Cd — Cadmium</option>
    <option value="V">V — Vanadium</option>
  </optgroup>
  <optgroup label="Other">
    <option value="Al">Al — Aluminium</option>
    <option value="Pb">Pb — Lead</option>
  </optgroup>
  <optgroup label="REE">
    <option value="La">La — Lanthanum</option>
    <option value="Ce">Ce — Cerium</option>
  </optgroup>
  <option value="other">Other (set mol. weight manually)</option>`;

// Map element symbols/names to their selector values (longest match wins)
const ELEM_DETECT = [
  {keys: ['nickel','Ni'],               val: 'Ni'},
  {keys: ['cobalt','Co'],               val: 'Co'},
  {keys: ['copper','Cu'],               val: 'Cu'},
  {keys: ['vanadium','Va','\bV\b'],   val: 'V'},
  {keys: ['zinc','Zn'],                 val: 'Zn'},
  {keys: ['cadmium','Cd'],              val: 'Cd'},
  {keys: ['manganese','Mn'],            val: 'Mn'},
  {keys: ['iron','Fe'],                 val: 'Fe'},
  {keys: ['alumin','Al'],               val: 'Al'},
  {keys: ['lead','Pb'],                 val: 'Pb'},
  {keys: ['lanthanum','La'],            val: 'La'},
  {keys: ['cerium','Ce'],               val: 'Ce'},
];

function detectElemFromCol(colName) {
  if (!colName) return null;
  const c = colName.trim();
  for (const e of ELEM_DETECT) {
    for (const k of e.keys) {
      // exact case-insensitive match OR word-boundary match
      if (c.toLowerCase() === k.toLowerCase()) return e.val;
      const re = new RegExp('(^|[^a-zA-Z])' + k + '([^a-zA-Z]|$)', 'i');
      if (re.test(c)) return e.val;
    }
  }
  return null;
}

function applyElemDetection(idx) {
  const col = teRowData[idx]?.col;
  if (!col) return;
  const elem = detectElemFromCol(col);
  if (!elem) return;
  const prefix = `te${idx + 1}`;
  const sel = document.getElementById(`${prefix}_elem`);
  if (sel && sel.value !== elem) {
    sel.value = elem;
    updateMolWt(prefix, elem);
  }
}

function teDepthChanged(val) {
  teDepthCol = val;
  teRawDepth = (teRawData[val] || []).map(Number).filter(isFinite);
  detectDepthUnit(val);
  renderTERows();
  updatePreviewChart();
}

function teDepthUnitChanged(val) {
  teDepthUnit = val;
  updatePreviewChart();
  renderAgePlot();
}

function detectDepthUnit(colName) {
  // 1. Check column header name for explicit unit
  const cl = (colName || '').toLowerCase();
  if (cl.includes('(mm)') || cl.includes('_mm') || cl.endsWith('mm') || cl.includes('mm_') || cl.includes(' mm')) {
    teDepthUnit = 'mm';
  } else if (cl.includes('(cm)') || cl.includes('_cm') || cl.endsWith('cm') || cl.includes('cm_') || cl.includes(' cm')) {
    teDepthUnit = 'cm';
  } else {
    // 2. Heuristic: max depth > 500 likely mm
    const maxDep = teRawDepth.length ? teRawDepth.reduce((a,b) => a > b ? a : b, 0) : 0;
    teDepthUnit = maxDep > 500 ? 'mm' : 'cm';
  }
  // Sync the dropdown
  const sel = document.getElementById('te-depth-unit-sel');
  if (sel) sel.value = teDepthUnit;
}

function addTERow() {
  const nonDepth = teColumns.filter(c => c !== teDepthCol);
  const used = teRowData.map(r => r.col);
  const next = nonDepth.find(c => !used.includes(c)) || nonDepth[0] || '';
  teRowData.push({col: next, unit: 'ppm'});
  renderTERows();
  renderTEParamCards();
  applyElemDetection(teRowData.length - 1);
  updatePreviewChart();
}

function removeTERow(idx) {
  if (teRowData.length <= 1) return;
  teRowData.splice(idx, 1);
  renderTERows();
  renderTEParamCards();
  updatePreviewChart();
}

function teSyncRow(idx) {
  const rowEl = document.querySelector(`.te-row[data-idx="${idx}"]`);
  if (!rowEl) return;
  teRowData[idx].col  = rowEl.querySelector('.te-col-sel').value;
  teRowData[idx].unit = rowEl.querySelector('.te-unit-sel').value;
  applyElemDetection(idx);
  updatePreviewChart();
  populateScatterSelectors();
}

function renderTERows() {
  const wrap = document.getElementById('te-rows');
  if (!wrap) return;
  const nonDepth = teColumns.filter(c => c !== teDepthCol);
  const colOpts = c => nonDepth.map(x =>
    `<option value="${x}"${x===c?' selected':''} >${x}</option>`).join('');
  const unitOpts = u => UNIT_OPTS.map(o =>
    `<option value="${o.v}"${o.v===u?' selected':''} >${o.l}</option>`).join('');
  wrap.innerHTML = teRowData.map((r, i) => `
    <div class="te-row" data-idx="${i}"
         style="display:grid;grid-template-columns:1fr 160px 28px;gap:8px;align-items:end;margin-bottom:8px">
      <div>
        <label style="font-size:11px;color:var(--muted)">TE${i+1} proxy column</label>
        <select class="te-col-sel" data-idx="${i}" onchange="teSyncRow(${i})"
                style="width:100%;margin-top:4px">${colOpts(r.col)}</select>
      </div>
      <div>
        <label style="font-size:11px;color:var(--muted)">Units</label>
        <select class="te-unit-sel" data-idx="${i}" onchange="teSyncRow(${i})"
                style="width:100%;margin-top:4px">${unitOpts(r.unit)}</select>
      </div>
      <button onclick="removeTERow(${i})"
              style="background:none;border:1px solid var(--border);color:var(--muted);
                     border-radius:4px;cursor:pointer;font-size:14px;height:32px;width:28px;
                     line-height:1;padding:0;align-self:end"
              title="Remove" ${teRowData.length<=1?'disabled':''}>×</button>
    </div>`).join('');
  populateScatterSelectors();
}

// ── Element-specific default parameters ───────────────────────────────
// Kd values from Lindeman et al. GCA 2022 (cave-analogue, +NOM conditions).
// Kd_mn = ln(Kd_NOM); users refine using the "implied ln(Kd)" hint.
// aq_conc are order-of-magnitude starting points — cave-specific.
const ELEM_DEFAULTS = {
  'Ni': { Kd_mn: -3.540, Kd_sd: 1.385, F: 0.00001, InertF: 0.10, aq_conc: 4.34, K_e: 1.0 },
  'Co': { Kd_mn: -0.891, Kd_sd: 1.385, F: 0.00001, InertF: 0.40, aq_conc: 0.47, K_e: 1.0 },
  'Cu': { Kd_mn: -0.083, Kd_sd: 1.385, F: 0.00001, InertF: 0.10, aq_conc: 1.000, K_e: 1.0 },
};
const ELEM_DEFAULT_FALLBACK = { Kd_mn: -2.000, Kd_sd: 1.385, F: 0.00001, InertF: 0.10, aq_conc: 1.000, K_e: 1.0 };

// Which element goes in which slot by default (before auto-detection)
const TE_POSITION_ELEM = { 1: 'Ni', 2: 'Co' };
const TE_POSITION_FALLBACK = 'Ni';

// Kp tooltip HTML — shown on first TE card as reference
const KP_TOOLTIP_HTML = `
  <div class="tooltip-box">
    <div class="tb-title">Theoretical Kp — Wang &amp; Xu (2001)</div>
    Setting Kp = −1 calculates the partition coefficient from cave temperature
    using the lattice strain model:<br><br>
    <code>log Kp = [ a(ΔGn_M − ΔGn_Ca) + b(r_M − r_Ca) − (ΔGf_M − ΔGf_Ca) ] / (−2.303·R·T)</code>
    where <em>r</em> is ionic radius relative to Ca²⁺, <em>ΔG<sub>f</sub></em> and
    <em>ΔG<sub>n</sub></em> are free energies of formation and hydration,
    <em>a</em> = 0.968 and <em>b</em> = 75.168 kcal/mol/Å are empirical constants,
    and <em>T</em> is cave temperature in Kelvin.<br><br>
    <strong style="color:var(--accent)">Speleothem-specific inorganic K<sub>d</sub> values</strong>
    (Lindeman et al., <em>GCA</em> 317, 2022 — cave-analogue crystal growth experiments
    under controlled pCO₂, temperature and humidity):
    <table style="width:100%; margin:6px 0; border-collapse:collapse; font-size:10px;">
      <tr style="color:var(--muted); border-bottom:1px solid var(--border)">
        <th style="text-align:left; padding:3px 6px">Element</th>
        <th style="text-align:right; padding:3px 6px">K<sub>d</sub> (inorganic)</th>
        <th style="text-align:right; padding:3px 6px">K<sub>d</sub> (+ NOM)</th>
        <th style="text-align:left; padding:3px 6px">Notes</th>
      </tr>
      <tr>
        <td style="padding:3px 6px; color:var(--texthi)">Co</td>
        <td style="text-align:right; padding:3px 6px">4.4</td>
        <td style="text-align:right; padding:3px 6px">0.41</td>
        <td style="padding:3px 6px; color:var(--muted)">PCP effect; NOM strongly suppresses</td>
      </tr>
      <tr style="background:rgba(255,255,255,0.03)">
        <td style="padding:3px 6px; color:var(--texthi)">Ni</td>
        <td style="text-align:right; padding:3px 6px">1.1</td>
        <td style="text-align:right; padding:3px 6px">0.029</td>
        <td style="padding:3px 6px; color:var(--muted)">PCP-insensitive (K<sub>d</sub> ≈ 1); recommended</td>
      </tr>
      <tr>
        <td style="padding:3px 6px; color:var(--texthi)">Cu</td>
        <td style="text-align:right; padding:3px 6px">44</td>
        <td style="text-align:right; padding:3px 6px">0.92</td>
        <td style="padding:3px 6px; color:var(--muted)">Very high K<sub>d</sub>; use with caution</td>
      </tr>
    </table>
    These Kp values are used internally by the model when Kp = −1. To override,
    enter a positive value directly. Note that in the presence of NOM, apparent K<sub>d</sub>
    values are far below the inorganic values — the model accounts for this via the
    X<sub>L</sub>, X<sub>F</sub>, and X<sub>I</sub> fractions.
  </div>`;

function makeTEParamCard(n, detectedElem) {
  const p = `te${n}`;
  // Use detected element if provided, otherwise fall back to position default
  const elem = detectedElem || TE_POSITION_ELEM[n] || TE_POSITION_FALLBACK;
  const d = ELEM_DEFAULTS[elem] || ELEM_DEFAULT_FALLBACK;
  const mw = MOL_WT[elem] || 58.693;
  const kp = KP_DEFAULT[elem] !== undefined ? KP_DEFAULT[elem] : -1;
  const labile = Math.max(0, 1 - d.F - d.InertF);

  // Build element selector with correct default selected
  const elemOpts = ELEM_OPTS.replace(
    `value="${elem}"`, `value="${elem}" selected`);
  // Kp field: include tooltip on TE1, Kp_theo display on all
  const showTooltip = (n === 1);
  const kpField = showTooltip
    ? `<div class="field has-tip">
        <div class="tip-label">
          <label>Partition coefficient Kp</label>
          <span class="tip-icon">?</span>
        </div>
        ${KP_TOOLTIP_HTML}
        <div style="display:flex; gap:6px; align-items:center">
          <input type="number" id="${p}_Kp" value="${kp}" step="0.01" style="flex:1"
                 oninput="updateTheoKp('${p}')">
          <input type="text" id="${p}_Kp_theo" readonly placeholder="theo. Kp"
                 title="Theoretical Kp from Wang &amp; Xu (2001), calculated when Kp = −1"
                 style="flex:1; opacity:0.55; cursor:not-allowed; font-size:10px; text-align:center">
        </div>
      </div>`
    : `<div class="field"><label>Partition coefficient Kp</label>
        <div style="display:flex; gap:6px; align-items:center">
          <input type="number" id="${p}_Kp" value="${kp}" step="0.01" style="flex:1"
                 oninput="updateTheoKp('${p}')">
          <input type="text" id="${p}_Kp_theo" readonly placeholder="theo. Kp"
                 title="Theoretical Kp from Wang &amp; Xu (2001), calculated when Kp = −1"
                 style="flex:1; opacity:0.55; cursor:not-allowed; font-size:10px; text-align:center">
        </div>
      </div>`;

  return `<div class="card" id="te-param-card-${n}">
    <div class="card-title">⚗️ Trace Element ${n}</div>
    <div class="form-grid">
      <div class="field"><label>Element</label>
        <select id="${p}_elem" onchange="updateMolWt('${p}', this.value)">${elemOpts}</select></div>
      <div class="field"><label>Molecular weight (g/mol)</label>
        <input type="number" id="${p}_mol_wt" value="${mw}" step="0.001"></div>
      ${kpField}
      <div class="field" style="grid-column:1/-1">
        <div style="display:flex;gap:6px;margin-bottom:8px">
          <button id="${p}_kd_mode_dr" onclick="setKdMode('${p}','driprate')"
                  class="btn btn-primary"
                  style="flex:1;padding:7px 10px;font-size:10px;height:auto">
            💧 From drip rate
          </button>
          <button id="${p}_kd_mode_lit" onclick="setKdMode('${p}','manual')"
                  class="btn btn-ghost"
                  style="flex:1;padding:7px 10px;font-size:10px;height:auto">
            📖 Literature / manual
          </button>
        </div>
        <div id="${p}_dr_wrap" style="margin-bottom:8px">
          <div style="margin-top:4px">
            <label style="font-size:11px;color:var(--muted)">
              Calibration window — anchor to most recent
              <strong id="${p}_cal_pct_label">5</strong>% of deposition
            </label>
            <div style="display:flex;gap:8px;align-items:center;margin-top:2px">
              <span style="font-size:9px;color:var(--muted)">1%</span>
              <input type="range" id="${p}_cal_pct" min="1" max="100" value="5" step="1"
                     style="flex:1;accent-color:var(--accent)"
                     oninput="document.getElementById('${p}_cal_pct_label').textContent=this.value;updateParamHints('${p}')">
              <span style="font-size:9px;color:var(--muted)">100%</span>
            </div>
            <div id="${p}_cal_hint" style="font-size:9.5px;color:var(--muted);margin-top:3px;min-height:13px;line-height:1.4"></div>
          </div>
          <div id="${p}_dr_hint" style="font-size:9.5px;color:var(--muted);margin-top:3px;min-height:13px;line-height:1.4"></div>
        </div>
        <div style="display:grid;grid-template-columns:1fr 1fr;gap:10px">
          <div>
            <label style="font-size:11px;color:var(--muted)">Mean ln(K<sub>d</sub>)</label>
            <input type="number" id="${p}_Kd_mn" value="${d.Kd_mn}" step="0.001"
                   oninput="updateParamHints('${p}')"
                   readonly style="opacity:0.6;cursor:not-allowed">
            <div id="${p}_Kd_mn_hint" style="font-size:9.5px;color:var(--muted);margin-top:3px;min-height:13px;line-height:1.4"></div>
          </div>
          <div>
            <label style="font-size:11px;color:var(--muted)">Std dev ln(K<sub>d</sub>)</label>
            <input type="number" id="${p}_Kd_sd" value="${d.Kd_sd}" step="0.001"
                   oninput="updateParamHints('${p}')">
            <div id="${p}_Kd_sd_hint" style="font-size:9.5px;color:var(--muted);margin-top:3px;min-height:13px;line-height:1.4"></div>
          </div>
        </div>
      </div>
      <div class="field">
        <label title="Partitioning efficiency multiplier on K₀. At K_e=1 the full thermodynamic Kp applies. Values <1 reduce effective partitioning (growth rate effects, surface saturation). Essentially Kp_effective = Kp × K_e — useful for keeping Kp at the published value while tuning empirically.">
          K<sub>e</sub> — partition efficiency
          <span style="color:var(--muted);font-weight:400;font-size:9px;cursor:help" title="K_e scales K₀ = Kp·Y_s·(Xa/Ya)·K_e. At 1.0 the thermodynamic Kp operates at full efficiency. Reduce below 1 to account for non-equilibrium growth, surface effects, or co-precipitation competition."> (?)</span>
        </label>
        <input type="number" id="${p}_K_e" value="${d.K_e}" step="0.1" min="0.001"
               oninput="updateParamHints('${p}')"></div>
      <div class="field"><label>X<sub>F</sub> — Fast fraction</label>
        <input type="number" id="${p}_F" value="${d.F}" step="0.001" min="0" max="1"
               oninput="updateLabile('${p}');updateParamHints('${p}')"></div>
      <div class="field"><label>X<sub>I</sub> — Inert fraction</label>
        <input type="number" id="${p}_InertF" value="${d.InertF}" step="0.001" min="0" max="1"
               oninput="updateLabile('${p}');updateParamHints('${p}')"></div>
      <div class="field"><label>X<sub>L</sub> — Labile (auto)</label>
        <input type="number" id="${p}_labile" value="${labile.toFixed(4)}" step="0.001" readonly
               style="opacity:0.6;cursor:not-allowed"></div>
      <div class="field fullonly" style="grid-column:1/-1">
        <label style="font-size:11px;font-weight:600;color:var(--texthi);margin-bottom:6px">
          Aqueous [<span id="${p}_aq_elem_label">${elem}</span>] concentration</label>
        <div style="display:flex;gap:6px;margin-bottom:8px">
          <button id="${p}-aq-mode-manual" class="mode-btn active"
                  onclick="setAqMode('${p}','manual')">Manual</button>
          <button id="${p}-aq-mode-csv" class="mode-btn"
                  onclick="setAqMode('${p}','csv')">Monitoring CSV</button>
        </div>
        <!-- Manual path -->
        <div id="${p}-aq-manual-block">
          <div style="display:grid;grid-template-columns:1fr 1fr;gap:8px">
            <div>
              <label style="font-size:10px;color:var(--muted)">Mean</label>
              <input type="number" id="${p}_aq_conc" value="${d.aq_conc}" step="0.001"
                     oninput="updateParamHints('${p}');fitAqPriorFromManual('${p}')">
            </div>
            <div>
              <label style="font-size:10px;color:var(--muted)">Std dev (optional)</label>
              <input type="number" id="${p}_aq_conc_sd" ${d.aq_conc_sd ? 'value="'+d.aq_conc_sd+'"' : ''} placeholder="blank = fixed"
                     step="0.001" oninput="fitAqPriorFromManual('${p}')">
            </div>
          </div>
          <select id="${p}_aq_unit" style="margin-top:4px;width:100%"
                  onchange="updateParamHints('${p}');fitAqPriorFromManual('${p}')">
            <option value="ppb" selected>ppb (µg/L)</option><option value="ppm">ppm (mg/L)</option>
          </select>
        </div>
        <!-- CSV path -->
        <div id="${p}-aq-csv-block" style="display:none">
          <div class="dropzone" id="dz-${p}_aq"
               ondragover="dzDrag(event,this)" ondragleave="dzLeave(this)"
               ondrop="dzDrop(event,'${p}_aq',this)">
            <input type="file" accept=".csv"
                   onchange="dzFile(event,'${p}_aq',this.parentElement)">
            <div class="dz-label">Drop <span class="${p}-aq-elem-dyn">${elem}</span> monitoring CSV</div>
            <div class="dz-name" id="name-${p}_aq"></div>
          </div>
          <div id="${p}-aq-col-selector" style="display:none;margin-top:8px">
            <label style="font-size:10px;color:var(--muted)">Concentration column</label>
            <select id="${p}_aq_csv_col" onchange="fitAqPriorFromCsv('${p}')"
                    style="width:100%;margin-bottom:4px"></select>
            <label style="font-size:10px;color:var(--muted)">Unit in file</label>
            <select id="${p}_aq_csv_unit" onchange="fitAqPriorFromCsv('${p}')" style="width:100%">
              <option value="ppb" selected>ppb (µg/L)</option><option value="ppm">ppm (mg/L)</option>
            </select>
          </div>
          <div style="display:flex;align-items:center;justify-content:space-between;margin-top:10px;margin-bottom:4px">
            <span style="font-size:10px;color:var(--muted)">Distribution fit</span>
            <button id="${p}-dist-chart_toggle" onclick="toggleLogScale('${p}-dist-chart')"
                    style="background:none;border:1px solid var(--border);color:var(--muted);
                           border-radius:4px;cursor:pointer;font-size:10px;padding:2px 8px">Log x</button>
          </div>
          <div style="height:160px;display:none" id="${p}-dist-chart-wrap">
            <canvas id="${p}-dist-chart"></canvas>
          </div>
          <div id="${p}-dist-desc" style="display:none;margin-top:8px;padding:8px 10px;
               background:var(--bg);border:1px solid var(--border);border-radius:6px;font-size:10px"></div>
        </div>
        <!-- Prior summary -->
        <div id="${p}-aq-prior-summary" style="font-size:10px;color:var(--muted);margin-top:6px;
             padding:6px 8px;background:var(--bg);border:1px solid var(--border);
             border-radius:4px;display:none"></div>
        <div id="${p}_aq_hint" style="font-size:9.5px;color:var(--muted);margin-top:3px;min-height:13px;line-height:1.4"></div>
      </div>
      <div class="field fullonly" style="grid-column:1/-1">
        <div id="${p}_pred_obs" style="font-size:10px;padding:6px 8px;background:var(--bg);border:1px solid var(--border);border-radius:4px;min-height:0"></div></div>
      <!-- h(V) vs drip rate diagnostic chart -->
      <div class="field fullonly" style="grid-column:1/-1">
        <div style="display:flex;align-items:center;justify-content:space-between;margin-bottom:4px">
          <span style="font-size:10px;color:var(--muted)">Forward model h(V) vs data</span>
        </div>
        <div style="height:180px">
          <canvas id="${p}_hv_chart"></canvas>
        </div>
      </div>
    </div></div>`;
}

function renderTEParamCards() {
  const n = teRowData.length;
  const container = document.getElementById('te-all-param-cards');
  if (!container) return;
  container.innerHTML = '';
  for (let i = 0; i < n; i++) {
    // Detect element from column name in the inputs section
    const detectedElem = detectElemFromCol(teRowData[i]?.col) || null;
    container.insertAdjacentHTML('beforeend', makeTEParamCard(i + 1, detectedElem));
    // Apply defaults for the detected element
    const finalElem = detectedElem || TE_POSITION_ELEM[i + 1] || TE_POSITION_FALLBACK;
    // Set the dropdown to the detected element
    const sel = document.getElementById(`te${i+1}_elem`);
    if (sel) sel.value = finalElem;
    updateMolWt(`te${i+1}`, finalElem);
    updateLabile(`te${i+1}`);
    updateTheoKp(`te${i+1}`);
  }
  updateAllParamHints();
  // Re-apply analysis mode to show/hide .fullonly elements in new cards
  setAnalysisMode(analysisMode);
}

function populateColSelectors(key, columns) {
  if (key === 'trace_elem1') {
    // TE CSV: drive the new dynamic row UI
    teColumns = columns;
    // Populate depth selector
    const depSel = document.getElementById('te-depth-sel');
    if (depSel) {
      depSel.innerHTML = columns.map(c => `<option value="${c}">${c}</option>`).join('');
      // Auto-select a depth column
      const depGuess = ['depth','Depth','DEPTH','dist','Dist'];
      const hit = depGuess.find(h => columns.includes(h)) || columns[0];
      depSel.value = hit;
      teDepthCol = hit;
    }
    // Show the config section
    const cfg = document.getElementById('te-config');
    if (cfg) cfg.style.display = '';
    // Initialise with one row
    const nonDepth = columns.filter(c => c !== teDepthCol);
    teRowData = [{col: nonDepth[0] || columns[0], unit: 'ppm'}];
    renderTERows();
    renderTEParamCards();
    // Auto-detect element for each initial row after cards are in DOM
    teRowData.forEach((_, i) => applyElemDetection(i));
    updateAllParamHints();
    // Re-init preprocessing now that depth col is known
    if (teRawData && Object.keys(teRawData).length) {
      initPreprocessing(teRawData, teDepthCol, Object.keys(teRawData).length
        ? (teRawData[teDepthCol] || teRawData[Object.keys(teRawData)[0]] || []).length
        : teRawN);
    }
    return;
  }
  // Default: depth/age and isotope uploads
  const wrap = document.getElementById('cols-' + key);
  if (!wrap) return;
  wrap.classList.add('show');
  // Only populate selects whose id starts with "sel-" — not unit pickers
  wrap.querySelectorAll('select[id^="sel-"]').forEach(sel => {
    const role = sel.id.split('-').pop();
    sel.innerHTML = columns.map(c => `<option value="${c}">${c}</option>`).join('');

    // Auto-detect column: try exact match first, then case-insensitive substring
    const exactGuesses = {
      depth:   ['depth','Depth','DEPTH','dist','Dist','distance','Distance'],
      age:     ['age','Age','AGE','yr_bp','yBP','age_BP','Age_BP','cal_age','CalAge'],
      age_err: ['error','Error','err','Err','sigma','Sigma','age_err','uncertainty','sd','SD','1sd','2sd'],
      proxy:   columns.filter(c => !c.toLowerCase().includes('depth')),
    };
    const substringGuesses = {
      depth:   ['depth','dist'],
      age:     ['age','yr','bp','ka'],
      age_err: ['err','sigma','uncert','sd'],
    };

    let found = (exactGuesses[role] || []).find(h => columns.includes(h));
    // If no exact match, try case-insensitive substring (skip depth-like cols for age)
    if (!found && substringGuesses[role]) {
      const subs = substringGuesses[role];
      found = columns.find(c => {
        const cl = c.toLowerCase();
        if (role === 'age') return subs.some(s => cl.includes(s)) && !cl.includes('depth') && !cl.includes('dist');
        if (role === 'age_err') return subs.some(s => cl.includes(s)) && !cl.includes('depth');
        return subs.some(s => cl.includes(s));
      });
    }
    if (found) sel.value = found;
    storeCol(key, role, sel.value);
  });
}

// ── Age-depth plot ────────────────────────────────────────────────────────
let agePlotData = null;
let agePlotChart = null;

// ── Linear regression ────────────────────────────────────────────────────
function linReg(xs, ys) {
  const n = xs.length;
  if (n < 2) return null;
  const mx = xs.reduce((a,b)=>a+b)/n, my = ys.reduce((a,b)=>a+b)/n;
  let ss = 0, sp = 0;
  for (let i=0;i<n;i++){ss+=(xs[i]-mx)**2; sp+=(xs[i]-mx)*(ys[i]-my);}
  if (ss===0) return null;
  const slope=sp/ss, int=my-slope*mx;
  return {predict: x => slope*x + int, slopeLo: slope, slopeHi: slope};
}

// ── PCHIP monotone cubic spline (Fritsch-Carlson 1980) ───────────────────
// Guarantees monotonicity between knots — no oscillation through kinks.
// Extrapolates linearly beyond data range using end tangents.
function pchip(xs, ys) {
  const n = xs.length;
  if (n < 2) return null;
  if (n === 2) {
    const s = (ys[1]-ys[0])/(xs[1]-xs[0]);
    return {predict: x => ys[0] + s*(x-xs[0]), slopeLo: s, slopeHi: s};
  }

  // 1. Secant slopes between knots
  const h = [], delta = [];
  for (let i=0;i<n-1;i++) {
    h[i]     = xs[i+1]-xs[i];
    delta[i] = (ys[i+1]-ys[i])/h[i];
  }

  // 2. Fritsch-Carlson tangent slopes
  const m = new Array(n);
  // Endpoints: one-sided
  m[0]   = delta[0];
  m[n-1] = delta[n-2];
  // Interior: harmonic mean, set to 0 at sign changes
  for (let i=1;i<n-1;i++) {
    if (delta[i-1]*delta[i] <= 0) {
      m[i] = 0;
    } else {
      const w1 = 2*h[i]+h[i-1], w2 = h[i]+2*h[i-1];
      m[i] = (w1+w2) / (w1/delta[i-1] + w2/delta[i]);
    }
  }

  // 3. Enforce monotonicity (Fritsch-Carlson condition)
  for (let i=0;i<n-1;i++) {
    if (Math.abs(delta[i]) < 1e-12) { m[i]=m[i+1]=0; continue; }
    const alpha=m[i]/delta[i], beta=m[i+1]/delta[i];
    if (alpha**2+beta**2 > 9) {
      const tau=3/Math.sqrt(alpha**2+beta**2);
      m[i]=tau*alpha*delta[i]; m[i+1]=tau*beta*delta[i];
    }
  }

  // 4. Evaluate cubic Hermite at x
  function predict(x) {
    if (x <= xs[0])   return ys[0]   + m[0]   * (x - xs[0]);
    if (x >= xs[n-1]) return ys[n-1] + m[n-1] * (x - xs[n-1]);
    // Binary search for segment
    let lo=0, hi=n-2;
    while (lo<hi) { const mid=(lo+hi)>>1; xs[mid+1]<=x ? lo=mid+1 : hi=mid; }
    const t=(x-xs[lo])/h[lo];
    const t2=t*t, t3=t2*t;
    return (2*t3-3*t2+1)*ys[lo] + (t3-2*t2+t)*h[lo]*m[lo]
         + (-2*t3+3*t2)*ys[lo+1] + (t3-t2)*h[lo]*m[lo+1];
  }

  return {predict, slopeLo: m[0], slopeHi: m[n-1]};
}

function renderAgePlot() {
  if (!agePlotData) return;
  const depCol = document.getElementById('sel-depth_age-depth')?.value;
  const ageCol = document.getElementById('sel-depth_age-age')?.value;
  const errCol = document.getElementById('sel-depth_age-age_err')?.value;
  if (!depCol || !ageCol) return;

  const rawD = agePlotData[depCol] || [];
  const rawA = agePlotData[ageCol] || [];
  const rawE = errCol ? (agePlotData[errCol] || []) : [];

  // Sanity check: if age values look like depths (max < 200) and depth values
  // look like ages (max > 200), the columns are probably swapped
  const maxA = rawA.reduce((m,v) => { const n=parseFloat(v); return isFinite(n)&&n>m?n:m; }, 0);
  const maxD = rawD.reduce((m,v) => { const n=parseFloat(v); return isFinite(n)&&n>m?n:m; }, 0);
  let useD = rawD, useA = rawA;
  if (maxA < 200 && maxD > 200) {
    console.warn('Age-depth auto-detect: columns appear swapped (age max=' + maxA + ', depth max=' + maxD + '). Swapping.');
    useD = rawA; useA = rawD;
  }

  const pts = useD.map((d,i) => ({d: parseFloat(d), a: parseFloat(useA[i]),
                                   e: parseFloat(rawE[i]) || 0}))
                  .filter(p => isFinite(p.d) && isFinite(p.a) && p.d >= 0);
  if (pts.length === 0) return;

  pts.sort((a,b) => a.d - b.d);
  const depths = pts.map(p=>p.d), ages = pts.map(p=>p.a);

  // Fit age=f(depth), then display as (age, depth) with depth on Y inverted
  const fitType = document.getElementById('age-fit-type')?.value || 'pchip';
  const fit = fitType === 'linear' ? linReg(depths, ages) : pchip(depths, ages);

  const dMin = Math.max(0, depths[0]), dMax = depths[depths.length-1];
  const dSpan = dMax - dMin || 1;
  const dExtMax = dMax + dSpan * 0.05;
  const N_CURVE = 300;
  const curvePts = [];
  for (let i = 0; i < N_CURVE; i++) {
    const d = dMin + i * (dExtMax - dMin) / (N_CURVE - 1);
    const a = fit ? fit.predict(d) : null;
    if (a !== null && isFinite(a)) curvePts.push({x: a, y: d});
  }

  const ageDataMin = ages.reduce((a,b) => a < b ? a : b, Infinity);
  const ageDataMax = ages.reduce((a,b) => a > b ? a : b, -Infinity);

  let ageMin = ageDataMin, ageMax = ageDataMax;
  if (fit) {
    const predMax = fit.predict(dExtMax);
    if (isFinite(predMax)) ageMax = Math.round(Math.max(ageMax, predMax));
    ageMin = Math.round(ageDataMin);
  }

  const extMin = document.getElementById('age-extrap-min');
  const extMax = document.getElementById('age-extrap-max');
  if (extMin && !extMin.dataset.manual) extMin.value = ageMin;
  if (extMax && !extMax.dataset.manual) extMax.value = ageMax;
  applyAgeRange(ageMin, ageMax);

  const wrap = document.getElementById('age-plot-wrap');
  if (wrap) wrap.style.display = '';
  const ctx = document.getElementById('agePlotCanvas').getContext('2d');
  if (agePlotChart) agePlotChart.destroy();

  const datasets = [];

  // 2σ error envelope — polygon approach
  if (errCol && agePlotData[errCol] && pts.some(p => p.e > 0)) {
    const ePts = pts.filter(p => p.e > 0).sort((a,b) => a.d - b.d);
    const envPts = [];
    ePts.forEach(p => envPts.push({x: p.a + p.e * 2, y: p.d}));
    [...ePts].reverse().forEach(p => envPts.push({x: p.a - p.e * 2, y: p.d}));
    if (envPts.length >= 3) {
      envPts.push({...envPts[0]}); // close polygon
      datasets.push({
        label: '2σ envelope',
        data: envPts, type: 'line',
        borderColor: 'rgba(76,201,160,0.2)',
        backgroundColor: 'rgba(76,201,160,0.08)',
        borderWidth: 1, pointRadius: 0, fill: true, tension: 0.2,
      });
    }
  }


  // Fit line (solid within data, dashed extrapolation)
  if (curvePts.length > 1) {
    const solidPts = curvePts.filter(p => p.y >= dMin && p.y <= dMax);
    const extPts   = curvePts.filter(p => p.y > dMax);
    if (solidPts.length && extPts.length) extPts.unshift(solidPts[solidPts.length-1]);

    datasets.push({
      label: fitType === 'linear' ? 'Linear fit' : 'PCHIP spline',
      data: solidPts, type: 'line',
      borderColor: 'rgba(76,201,160,0.7)', borderWidth: 1.8,
      pointRadius: 0, fill: false, tension: 0,
    });
    if (extPts.length > 1) {
      datasets.push({
        label: 'Extrapolation', data: extPts, type: 'line',
        borderColor: 'rgba(180,180,180,0.5)', borderWidth: 1.5,
        borderDash: [4,3], pointRadius: 0, fill: false, tension: 0,
      });
    }
  }

  // Error bars (horizontal — age uncertainty at each depth)
  if (errCol && agePlotData[errCol]) {
    const errSegs = [];
    for (const p of pts) {
      if (p.e > 0) {
        errSegs.push({x: p.a - p.e, y: p.d}, {x: p.a + p.e, y: p.d}, {x: NaN, y: NaN});
      }
    }
    if (errSegs.length) {
      datasets.push({
        label: 'Age error (\u00b11\u03c3)',
        data: errSegs, type: 'line',
        borderColor: 'rgba(247,164,64,0.5)', borderWidth: 1.5,
        pointRadius: 0, fill: false, spanGaps: false,
      });
    }
  }

  // Dated points
  datasets.push({
    label: 'Dated points',
    data: pts.map(p => ({x: p.a, y: p.d})),
    type: 'scatter',
    backgroundColor: 'rgba(76,201,160,0.85)',
    pointRadius: 5, pointHoverRadius: 7,
  });

  // Age-range shading + hiatus exclusion zones (both on X axis)
  const shadingPlugin = {
    id: 'ageRangeShading',
    afterDraw(chart) {
      const {ctx: c, chartArea: {left,right,top,bottom}, scales: {x}} = chart;
      // Age range band
      const mn = parseFloat(document.getElementById('calage_min')?.value);
      const mx = parseFloat(document.getElementById('calage_max')?.value);
      if (isFinite(mn) && isFinite(mx)) {
        const xMn = x.getPixelForValue(mn), xMx = x.getPixelForValue(mx);
        const xL = Math.min(xMn, xMx), xR = Math.max(xMn, xMx);
        c.save();
        c.fillStyle = 'rgba(76,201,160,0.07)';
        c.fillRect(xL, top, xR-xL, bottom-top);
        c.strokeStyle = 'rgba(76,201,160,0.4)';
        c.setLineDash([3,3]); c.lineWidth = 1;
        c.beginPath(); c.moveTo(xMn,top); c.lineTo(xMn,bottom); c.stroke();
        c.beginPath(); c.moveTo(xMx,top); c.lineTo(xMx,bottom); c.stroke();
        c.restore();
      }
      // Hiatus exclusion zones
      if (typeof hiatusZones !== 'undefined' && hiatusZones.length > 0) {
        c.save();
        hiatusZones.forEach(z => {
          const xF = x.getPixelForValue(z.from), xT = x.getPixelForValue(z.to);
          const xL = Math.min(xF, xT), xR = Math.max(xF, xT);
          c.fillStyle = 'rgba(255,80,80,0.12)';
          c.fillRect(xL, top, xR-xL, bottom-top);
          c.strokeStyle = 'rgba(255,80,80,0.4)';
          c.setLineDash([2,2]); c.lineWidth = 1;
          c.beginPath(); c.moveTo(xL,top); c.lineTo(xL,bottom); c.stroke();
          c.beginPath(); c.moveTo(xR,top); c.lineTo(xR,bottom); c.stroke();
          // Label
          c.fillStyle = 'rgba(255,80,80,0.6)';
          c.font = '9px Arial, sans-serif';
          c.textAlign = 'center';
          c.fillText('hiatus', (xL+xR)/2, top + 12);
        });
        c.restore();
      }
    }
  };

  try {
    agePlotChart = new Chart(ctx, {
      data: {datasets},
      options: {
        animation: false, responsive: true, maintainAspectRatio: false,
        plugins: {
          legend: { labels: { color:'#8fa4b5', font:{size:10},
            filter: item => !(item.text || '').startsWith('_') } },
          tooltip: { callbacks: {
            label: item => {
              const y = item.parsed?.y, x = item.parsed?.x;
              if (!isFinite(y) || !isFinite(x)) return '';
              return y.toFixed(2) + ' ' + teDepthUnit + ' depth \u2192 ' + x.toFixed(0) + ' yrs BP';
            }
          }}
        },
        scales: {
          x: {
            type: 'linear',
            title: {display:true, text:'Age (yrs BP)', color:'#8fa4b5', font:{size:11}},
            grid: {color:'rgba(255,255,255,0.05)'},
            ticks: {color:'#8fa4b5', font:{size:10}},
          },
          y: {
            type: 'linear',
            reverse: true,
            min: 0,
            title: {display:true, text:'Depth (' + teDepthUnit + ')', color:'#8fa4b5', font:{size:11}},
            grid: {color:'rgba(255,255,255,0.05)'},
            ticks: {color:'#8fa4b5', font:{size:10}},
          }
        }
      },
      plugins: [shadingPlugin],
    });
  } catch(e) { console.error('Age plot chart error:', e); }
  // Update growth rate chart
  updateGrowthRate();
}


function applyAgeRange(minAge, maxAge) {
  const lo = Math.min(minAge, maxAge), hi = Math.max(minAge, maxAge);
  document.getElementById('calage_min').value = lo;
  document.getElementById('calage_max').value = hi;
}

function onExtrapChange() {
  const extMin = document.getElementById('age-extrap-min');
  const extMax = document.getElementById('age-extrap-max');
  if (extMin) extMin.dataset.manual = '1';
  if (extMax) extMax.dataset.manual = '1';
  const mn = parseFloat(extMin?.value), mx = parseFloat(extMax?.value);
  if (isFinite(mn) && isFinite(mx)) {
    applyAgeRange(mn, mx);
    if (agePlotChart) agePlotChart.update();
  }
}

function syncExtrapFromFields() {
  // When user edits calage_min/max directly, reflect back into extrap inputs
  const mn = document.getElementById('calage_min')?.value;
  const mx = document.getElementById('calage_max')?.value;
  const extMin = document.getElementById('age-extrap-min');
  const extMax = document.getElementById('age-extrap-max');
  if (extMin) { extMin.value = mn; extMin.dataset.manual = '1'; }
  if (extMax) { extMax.value = mx; extMax.dataset.manual = '1'; }
  if (agePlotChart) agePlotChart.update();
}

// ── Analysis mode ────────────────────────────────────────────────────────
let analysisMode = 'full';

// ── Growth rate & hiatus detection ──────────────────────────────────────
let growthRateChart = null;
let hiatusZones = [];  // [{from: ageMin, to: ageMax}, ...]
let _lastGrowthData = null;  // {ages: [], rates: []}

function computeGrowthRate() {
  // Compute growth rate (mm/yr) from the current age-depth fit
  if (!agePlotData) return null;
  const depCol = document.getElementById('sel-depth_age-depth')?.value;
  const ageCol = document.getElementById('sel-depth_age-age')?.value;
  if (!depCol || !ageCol) return null;

  const rawD = agePlotData[depCol] || [];
  const rawA = agePlotData[ageCol] || [];
  const maxA = rawA.reduce((m,v) => { const n=parseFloat(v); return isFinite(n)&&n>m?n:m; }, 0);
  const maxD = rawD.reduce((m,v) => { const n=parseFloat(v); return isFinite(n)&&n>m?n:m; }, 0);
  let useD = rawD, useA = rawA;
  if (maxA < 200 && maxD > 200) { useD = rawA; useA = rawD; }

  const pts = useD.map((d,i) => ({d: parseFloat(d), a: parseFloat(useA[i])}))
                  .filter(p => isFinite(p.d) && isFinite(p.a) && p.d >= 0);
  if (pts.length < 2) return null;
  pts.sort((a,b) => a.d - b.d);
  const depths = pts.map(p=>p.d), ages = pts.map(p=>p.a);

  const fitType = document.getElementById('age-fit-type')?.value || 'pchip';
  const fit = fitType === 'linear' ? linReg(depths, ages) : pchip(depths, ages);
  if (!fit) return null;

  // Compute growth rate at regular age intervals
  const ageMin = ages.reduce((a,b)=>a<b?a:b, Infinity);
  const ageMax = ages.reduce((a,b)=>a>b?a:b, -Infinity);
  const nPts = 200;
  const dMin = Math.max(0, depths[0]), dMax = depths[depths.length-1];
  const result = {ages: [], rates: [], depths: []};

  for (let i = 0; i < nPts; i++) {
    const d = dMin + i * (dMax - dMin) / (nPts - 1);
    const age = fit.predict(d);
    // Numerical derivative: dDepth/dAge ≈ Δd/Δa
    const dd = 0.01;
    const ageP = fit.predict(d + dd);
    const dAge = ageP - age;
    // Growth rate in mm/yr (depth unit assumed cm unless mm detected)
    const depthMM = teDepthUnit === 'mm' ? dd : dd * 10;
    const rateMMperYr = dAge !== 0 ? Math.abs(depthMM / dAge) : 0;
    result.ages.push(age);
    result.rates.push(rateMMperYr);
    result.depths.push(d);
  }
  return result;
}

function updateGrowthRate() {
  _lastGrowthData = computeGrowthRate();
  renderGrowthRateChart();
}

function renderGrowthRateChart() {
  const data = _lastGrowthData;
  const canvas = document.getElementById('growthRateCanvas');
  if (!canvas || !data) return;
  if (growthRateChart) growthRateChart.destroy();

  const threshold = parseFloat(document.getElementById('hiatus-threshold')?.value) || 0.005;
  const ratePts = data.ages.map((a, i) => ({x: a, y: data.rates[i]}));

  // Highlight below-threshold regions
  const belowPts = ratePts.map(p => ({x: p.x, y: p.y < threshold ? p.y : null}));

  const datasets = [
    { label: 'Growth rate',
      data: ratePts, type: 'line',
      borderColor: 'rgba(76,201,160,0.8)', borderWidth: 1.5,
      pointRadius: 0, fill: false, tension: 0.3, yAxisID: 'y' },
    { label: 'Below threshold',
      data: belowPts, type: 'line',
      borderColor: 'rgba(255,80,80,0.7)', backgroundColor: 'rgba(255,80,80,0.15)',
      borderWidth: 1.5, pointRadius: 0, fill: true, tension: 0.3,
      spanGaps: false, yAxisID: 'y' },
    { label: 'Threshold',
      data: [{x: data.ages[0], y: threshold}, {x: data.ages[data.ages.length-1], y: threshold}],
      type: 'line', borderColor: 'rgba(255,160,50,0.5)', borderWidth: 1,
      borderDash: [4,3], pointRadius: 0, fill: false, yAxisID: 'y' },
  ];

  // Highlight exclusion zones
  hiatusZones.forEach((z, zi) => {
    datasets.push({
      label: zi === 0 ? 'Exclusion zones' : '_ez' + zi,
      data: [{x: z.from, y: 0}, {x: z.from, y: threshold * 3}, {x: z.to, y: threshold * 3}, {x: z.to, y: 0}],
      type: 'line', borderColor: 'rgba(255,80,80,0.3)',
      backgroundColor: 'rgba(255,80,80,0.08)',
      borderWidth: 1, pointRadius: 0, fill: true, tension: 0,
    });
  });

  growthRateChart = new Chart(canvas, {
    data: {datasets},
    options: {
      animation: false, responsive: true, maintainAspectRatio: false,
      plugins: { legend: { labels: { color:'#8fa4b5', font:{size:10},
        filter: item => !(item.text || '').startsWith('_') } } },
      scales: {
        x: { type: 'linear',
             title: {display:true, text:'Age (yrs BP)', color:'#8fa4b5', font:{size:10}},
             ticks: {color:'#8fa4b5', font:{size:10}},
             grid: {color:'rgba(255,255,255,0.04)'} },
        y: { type: 'linear', min: 0,
             title: {display:true, text:'Growth rate (mm/yr)', color:'#8fa4b5', font:{size:10}},
             ticks: {color:'#8fa4b5', font:{size:10}},
             grid: {color:'rgba(255,255,255,0.04)'} },
      }
    }
  });
}

function autoDetectHiatuses() {
  const data = _lastGrowthData || computeGrowthRate();
  if (!data) return;
  const threshold = parseFloat(document.getElementById('hiatus-threshold')?.value) || 0.005;

  // Find contiguous age intervals where growth rate < threshold
  hiatusZones = [];
  let inHiatus = false, startAge = 0;
  for (let i = 0; i < data.ages.length; i++) {
    const below = data.rates[i] < threshold;
    if (below && !inHiatus) {
      inHiatus = true;
      startAge = data.ages[i];
    } else if (!below && inHiatus) {
      inHiatus = false;
      hiatusZones.push({from: Math.round(startAge), to: Math.round(data.ages[i])});
    }
  }
  if (inHiatus) {
    hiatusZones.push({from: Math.round(startAge), to: Math.round(data.ages[data.ages.length-1])});
  }

  renderHiatusZoneList();
  renderGrowthRateChart();
  // Also re-render the age-depth chart to show zone overlays
  renderAgePlot();
}

function addHiatusZone() {
  hiatusZones.push({from: 0, to: 1000});
  renderHiatusZoneList();
  renderGrowthRateChart();
}

function removeHiatusZone(idx) {
  hiatusZones.splice(idx, 1);
  renderHiatusZoneList();
  renderGrowthRateChart();
  renderAgePlot();
}

function updateHiatusZone(idx, field, value) {
  hiatusZones[idx][field] = parseFloat(value) || 0;
  renderGrowthRateChart();
  renderAgePlot();
}

function renderHiatusZoneList() {
  const list = document.getElementById('hiatus-zone-list');
  const hint = document.getElementById('hiatus-zone-hint');
  if (!list) return;

  if (hiatusZones.length === 0) {
    list.innerHTML = '';
    if (hint) hint.style.display = '';
    return;
  }
  if (hint) hint.style.display = 'none';

  list.innerHTML = hiatusZones.map((z, i) => `
    <div style="display:grid;grid-template-columns:1fr 1fr 28px;gap:6px;align-items:end;margin-bottom:6px">
      <div>
        <label style="font-size:10px;color:var(--muted)">From (yrs BP)</label>
        <input type="number" value="${z.from}" style="width:100%"
               oninput="updateHiatusZone(${i},'from',this.value)">
      </div>
      <div>
        <label style="font-size:10px;color:var(--muted)">To (yrs BP)</label>
        <input type="number" value="${z.to}" style="width:100%"
               oninput="updateHiatusZone(${i},'to',this.value)">
      </div>
      <button onclick="removeHiatusZone(${i})"
              style="background:none;border:1px solid var(--border);color:var(--muted);
                     border-radius:4px;cursor:pointer;font-size:14px;height:32px;width:28px;
                     line-height:1;padding:0;align-self:end" title="Remove">×</button>
    </div>`).join('');
}

function setAnalysisMode(mode) {
  analysisMode = mode;
  const isFull = mode === 'full';
  // Toggle button styles
  document.getElementById('mode-btn-full').className = isFull ? 'btn btn-primary' : 'btn btn-ghost';
  document.getElementById('mode-btn-semi').className = isFull ? 'btn btn-ghost' : 'btn btn-primary';
  // Callout text
  document.getElementById('mode-callout-full').style.display = isFull ? '' : 'none';
  document.getElementById('mode-callout-semi').style.display = isFull ? 'none' : '';
  // Show/hide all .fullonly elements
  document.querySelectorAll('.fullonly').forEach(el => {
    el.style.display = isFull ? '' : 'none';
  });
  // Update Smart & Friedrich tab availability
  const sfTab = document.getElementById('tab-sf');
  if (sfTab) {
    sfTab.style.opacity = isFull ? '1' : '0.4';
    sfTab.title = isFull ? '' : 'Smart & Friedrich classification requires Full Quantification mode';
  }
}

// ── TE preprocessing state ───────────────────────────────────────────────
let teRawData    = {};   // {colName: [values...]} — raw uploaded data
let teRawDepth   = [];   // raw depth column values
let teRawN       = 0;    // original row count
let tePreviewChart = null;

function initPreprocessing(rawData, depthCol, rowCount) {
  teRawData  = rawData;
  teRawDepth = (rawData[depthCol] || []).map(Number).filter(isFinite);
  teRawN     = rowCount;
  detectDepthUnit(depthCol);
  document.getElementById('te-target-n').value  = Math.round(rowCount / 5);
  document.getElementById('te-target-n').max    = rowCount;
  document.getElementById('te-row-badge').textContent = `${rowCount} rows uploaded`;
  document.getElementById('te-preproc').style.display = '';
  document.getElementById('te-preproc-status').textContent = 'No preprocessing applied — using raw data.';
  updatePreviewChart();
  updateRuntimeEstimate();
}

// Windowed sigma-clip: for each point compute local mean+SD from a centred
// window of `win` neighbours (0 = global). Returns bool[] of outlier flags.
function windowedSigmaClip(arr, sigma, win) {
  const n = arr.length;
  const flags = new Array(n).fill(false);
  if (win <= 0) {
    // Global clip
    const valid = arr.filter(isFinite);
    if (!valid.length) return flags;
    const mu = valid.reduce((a,b)=>a+b,0)/valid.length;
    const sd = Math.sqrt(valid.map(v=>(v-mu)**2).reduce((a,b)=>a+b,0)/valid.length);
    return arr.map(v => isFinite(v) && Math.abs(v - mu) > sigma * sd);
  }
  const half = Math.floor(win / 2);
  for (let i = 0; i < n; i++) {
    if (!isFinite(arr[i])) continue;
    const lo = Math.max(0, i - half), hi = Math.min(n, i + half + 1);
    const slice = arr.slice(lo, hi).filter(isFinite);
    if (slice.length < 3) continue;
    const mu = slice.reduce((a,b)=>a+b,0)/slice.length;
    const sd = Math.sqrt(slice.map(v=>(v-mu)**2).reduce((a,b)=>a+b,0)/slice.length);
    if (sd > 0 && Math.abs(arr[i] - mu) > sigma * sd) flags[i] = true;
  }
  return flags;
}

const TE_CHART_COLORS = [
  {bg:'rgba(76,201,160,0.5)', border:'rgba(76,201,160,0.8)', outlier:'rgba(255,100,100,0.7)'},
  {bg:'rgba(247,164,64,0.5)', border:'rgba(247,164,64,0.8)', outlier:'rgba(200,80,40,0.7)'},
  {bg:'rgba(168,85,247,0.5)', border:'rgba(168,85,247,0.8)', outlier:'rgba(220,60,180,0.7)'},
  {bg:'rgba(96,165,250,0.5)', border:'rgba(96,165,250,0.8)', outlier:'rgba(60,100,200,0.7)'},
];

function updatePreviewChart() {
  // Show ALL configured TE columns stacked with different colors
  const dep = teRawDepth;
  if (!dep.length) return;

  const sigma   = parseFloat(document.getElementById('te-sigma').value) || 3;
  const targetN = parseInt(document.getElementById('te-target-n').value) || teRawN;
  const winSize = parseInt(document.getElementById('te-winsize').value) || 0;
  const step    = Math.max(1, Math.round(teRawN / targetN));

  const datasets = [];
  let totalOut = 0;
  const perTeStats = [];

  teRowData.forEach((r, ti) => {
    const col = r.col;
    if (!col || !teRawData[col]) return;
    const raw = teRawData[col].map(Number);
    if (raw.length !== dep.length) return;
    const isOutlier = windowedSigmaClip(raw, sigma, winSize);
    const c = TE_CHART_COLORS[ti % TE_CHART_COLORS.length];
    const elem = detectElemFromCol(col) || `TE${ti+1}`;
    const ktd = [], out = [];
    raw.forEach((v, i) => {
      if (!isFinite(v)) return;
      if (isOutlier[i]) { out.push({x: dep[i], y: v}); return; }
      if (i % step === 0) ktd.push({x: dep[i], y: v});
    });
    totalOut += out.length;
    perTeStats.push({elem, nOut: out.length, nKept: ktd.length});
    datasets.push({ label: `${elem} kept`, data: ktd,
      backgroundColor: c.bg, pointRadius: 2, pointHoverRadius: 4 });
    if (out.length) datasets.push({ label: `${elem} outlier`, data: out,
      backgroundColor: c.outlier, borderColor: c.outlier, pointRadius: 3,
      pointStyle: 'crossRot' });
  });

  if (!datasets.length) return;
  const ctx = document.getElementById('te-preview-chart').getContext('2d');
  if (tePreviewChart) tePreviewChart.destroy();
  tePreviewChart = new Chart(ctx, {
    type: 'scatter', data: { datasets },
    options: {
      animation: false, responsive: true, maintainAspectRatio: false,
      plugins: { legend: { labels: { color:'#8fa4b5', font:{size:10} } } },
      scales: {
        x: { ticks:{color:'#8fa4b5',font:{size:10}}, grid:{color:'rgba(255,255,255,0.04)'},
             title:{display:true,text:`Depth (${teDepthUnit})`,color:'#8fa4b5',font:{size:10}} },
        y: { ticks:{color:'#8fa4b5',font:{size:10}}, grid:{color:'rgba(255,255,255,0.04)'},
             title:{display:true,text:'Concentration',color:'#8fa4b5',font:{size:10}} }
      }
    }
  });
  const nFinal = Math.round((teRawN - totalOut) / Math.max(1, step));
  const perDetail = perTeStats.map(s => `${s.elem}: ${s.nOut} outliers`).join(' · ');
  document.getElementById('te-row-badge').innerHTML =
    `${teRawN} rows \u2192 ~${nFinal} after preprocessing<br>` +
    `<span style="font-size:9px">${perDetail}</span>`;
  updateRuntimeEstimate();
}

// ── TE scatter plot ──────────────────────────────────────────────────────
let teScatterChart = null;

function populateScatterSelectors() {
  const panel = document.getElementById('te-scatter-panel');
  if (!panel) return;
  if (teRowData.length < 2) { panel.style.display = 'none'; return; }
  panel.style.display = '';
  const selX = document.getElementById('te-scatter-x');
  const selY = document.getElementById('te-scatter-y');
  if (!selX || !selY) return;
  const opts = teRowData.map((r,i) => {
    const elem = detectElemFromCol(r.col) || `TE${i+1}`;
    return `<option value="${r.col}">${elem} (${r.col})</option>`;
  }).join('');
  selX.innerHTML = opts;
  selY.innerHTML = opts;
  if (teRowData.length >= 2) { selX.selectedIndex = 0; selY.selectedIndex = 1; }
  updateScatterPlot();
}

function updateScatterPlot() {
  const colX = document.getElementById('te-scatter-x')?.value;
  const colY = document.getElementById('te-scatter-y')?.value;
  const fitType = document.getElementById('te-scatter-fit')?.value || 'none';
  if (!colX || !colY || colX === colY || !teRawData[colX] || !teRawData[colY]) return;

  const xRaw = teRawData[colX].map(Number);
  const yRaw = teRawData[colY].map(Number);
  const pts = [];
  for (let i = 0; i < Math.min(xRaw.length, yRaw.length); i++) {
    if (isFinite(xRaw[i]) && isFinite(yRaw[i])) pts.push({x: xRaw[i], y: yRaw[i]});
  }
  if (pts.length < 3) return;

  // Compute Pearson r
  const xs = pts.map(p=>p.x), ys = pts.map(p=>p.y);
  const n = xs.length;
  const mx = xs.reduce((a,b)=>a+b,0)/n, my = ys.reduce((a,b)=>a+b,0)/n;
  let sxx=0, syy=0, sxy=0;
  for (let i=0;i<n;i++){sxx+=(xs[i]-mx)**2;syy+=(ys[i]-my)**2;sxy+=(xs[i]-mx)*(ys[i]-my);}
  const r = sxy / Math.sqrt(sxx*syy);
  const r2 = r*r;

  const datasets = [
    { label: 'Data', data: pts, type: 'scatter',
      backgroundColor: 'rgba(76,201,160,0.4)', pointRadius: 2, pointHoverRadius: 4 },
  ];

  // Fit line
  if (fitType === 'linear' && sxx > 0) {
    const slope = sxy/sxx, inter = my - slope*mx;
    const xMin = Math.min(...xs), xMax = Math.max(...xs);
    datasets.push({ label: `Linear (r\u00b2=${r2.toFixed(3)})`,
      data: [{x:xMin,y:slope*xMin+inter},{x:xMax,y:slope*xMax+inter}],
      type: 'line', borderColor:'#f7a440', borderWidth:2, pointRadius:0, fill:false });
  } else if (fitType === 'exp') {
    // ln(y) = a + b*x → y = exp(a)*exp(b*x)
    const lnys = ys.map(y => y > 0 ? Math.log(y) : NaN).filter(isFinite);
    const validPts = pts.filter(p => p.y > 0);
    if (validPts.length > 3) {
      const vx = validPts.map(p=>p.x), vly = validPts.map(p=>Math.log(p.y));
      const vmx = vx.reduce((a,b)=>a+b,0)/vx.length, vmy = vly.reduce((a,b)=>a+b,0)/vly.length;
      let vss=0,vsp=0;
      for(let i=0;i<vx.length;i++){vss+=(vx[i]-vmx)**2;vsp+=(vx[i]-vmx)*(vly[i]-vmy);}
      if (vss>0) {
        const b=vsp/vss, a=vmy-b*vmx;
        const xMin=Math.min(...vx), xMax=Math.max(...vx);
        const curve = Array.from({length:50},(_,i)=>{
          const x=xMin+i*(xMax-xMin)/49; return {x, y:Math.exp(a+b*x)};});
        datasets.push({ label: `Exp fit (r\u00b2=${r2.toFixed(3)})`,
          data: curve, type: 'line', borderColor:'#a855f7', borderWidth:2,
          pointRadius:0, fill:false, tension:0.3 });
      }
    }
  }

  const canvas = document.getElementById('te-scatter-chart');
  if (!canvas) return;
  if (teScatterChart) teScatterChart.destroy();
  const elemX = detectElemFromCol(colX) || colX;
  const elemY = detectElemFromCol(colY) || colY;
  teScatterChart = new Chart(canvas, {
    data: { datasets },
    options: {
      animation: false, responsive: true, maintainAspectRatio: false,
      plugins: { legend: { labels: { color:'#8fa4b5', font:{size:10} } } },
      scales: {
        x: { title:{display:true,text:elemX,color:'#8fa4b5',font:{size:11}},
             ticks:{color:'#8fa4b5',font:{size:10}}, grid:{color:'rgba(255,255,255,0.04)'} },
        y: { title:{display:true,text:elemY,color:'#8fa4b5',font:{size:11}},
             ticks:{color:'#8fa4b5',font:{size:10}}, grid:{color:'rgba(255,255,255,0.04)'} },
      }
    }
  });
  document.getElementById('te-scatter-stats').innerHTML =
    `n=${n} · Pearson r = <strong>${r.toFixed(3)}</strong> · r\u00b2 = <strong>${r2.toFixed(3)}</strong>` +
    (r > 0.5 ? ' · <span style="color:var(--accent)">positive correlation — consistent with shared drip rate control</span>'
     : r < -0.3 ? ' · <span style="color:#ffa032">negative correlation — different processes likely dominate</span>'
     : ' · <span style="color:var(--muted)">weak correlation — limited evidence for shared drip rate control</span>');
}

function applyPreprocessing() {
  const targetN = parseInt(document.getElementById('te-target-n').value) || teRawN;
  const sigma   = parseFloat(document.getElementById('te-sigma').value) || 3;
  const depCol  = teDepthCol;
  const st = document.getElementById('te-preproc-status');
  st.textContent = 'Applying…';
  fetch('/preprocess', {
    method: 'POST',
    headers: {'Content-Type':'application/json'},
    body: JSON.stringify({dep_col: depCol, target_n: targetN, sigma: sigma,
                           win_size: parseInt(document.getElementById('te-winsize').value) || 0})
  })
  .then(r => r.json())
  .then(d => {
    if (d.error) { st.textContent = 'Error: ' + d.error; return; }
    teRawN = d.rows;
    st.style.color = 'var(--accent)';
    st.textContent = `✓ Saved: ${d.rows} rows (was ${d.original_rows}; ${d.removed_outliers} outliers replaced)`;
    document.getElementById('te-row-badge').textContent = `${d.rows} rows (preprocessed)`;
    updateRuntimeEstimate();
  });
}

// ── Runtime estimate ─────────────────────────────────────────────────────
// Empirical: ~0.025 s per (depth point × calage step) on a typical laptop.
// Empirical constant: ~900s BayProX for 387pts × 9500yr × 2 passes (depth/age + 1 TE PDist).
// = 900 / (387 × 9500 × 2) ≈ 0.000122 s per point per year per pass.
// Hardware varies 2–5×; treat as order-of-magnitude guide.
const BAYPROX_SEC_PER_PT_YR = 0.000122;

function updateRuntimeEstimate() {
  const el = document.getElementById('runtime-estimate-body');
  if (!el) return;
  const nPts  = teRawN || 0;
  const cMin  = parseFloat(document.getElementById('calage_min')?.value) || 0;
  const cMax  = parseFloat(document.getElementById('calage_max')?.value) || 10000;
  const nAge  = Math.abs(cMax - cMin);
  const nTEs  = Math.max(teRowData.length, 1);
  const useCached = document.getElementById('use_cached_proxy')?.checked;

  if (nPts === 0) { el.innerHTML = '<span style="color:var(--muted)">Upload TE data to see estimate.</span>'; return; }

  // depth/age BayProX pass + one PDist pass per TE
  const nPasses    = 1 + nTEs;
  const bayproxSec = useCached ? 0 : nPts * nAge * BAYPROX_SEC_PER_PT_YR * nPasses;
  const dripSec    = nPts * nTEs * 0.05;  // fast
  const totalSec   = bayproxSec + dripSec + 60;  // +60s overhead

  const fmt = s => s < 120 ? `~${Math.round(s)}s`
                  : s < 3600 ? `~${Math.round(s/60)} min`
                  : `~${(s/3600).toFixed(1)} hr`;

  const warn = totalSec > 1200;
  const tip  = nPts > 500 && !useCached
    ? `<span style="color:var(--accent2)">💡 ${nPts} depth points is the main cost driver.
       Use the preprocessing panel to downsample, or check
       <em>Use cached proxy record</em> if you already ran BayProX.</span>`
    : '';

  el.innerHTML = `
    <div style="display:grid;grid-template-columns:repeat(3,1fr);gap:8px;margin-bottom:8px">
      <div style="background:var(--bg);border:1px solid var(--border);border-radius:6px;padding:8px 10px">
        <div style="font-size:10px;color:var(--muted)">Depth points</div>
        <div style="font-size:16px;color:var(--texthi);font-weight:600">${nPts.toLocaleString()}</div>
      </div>
      <div style="background:var(--bg);border:1px solid var(--border);border-radius:6px;padding:8px 10px">
        <div style="font-size:10px;color:var(--muted)">Age range</div>
        <div style="font-size:16px;color:var(--texthi);font-weight:600">${nAge.toLocaleString()} yr</div>
      </div>
      <div style="background:var(--bg);border:1px solid var(--border);border-radius:6px;padding:8px 10px;
                  ${warn ? 'border-color:rgba(255,160,50,0.5)' : ''}">
        <div style="font-size:10px;color:var(--muted)">Est. runtime</div>
        <div style="font-size:16px;font-weight:600;color:${warn ? '#ffa032' : 'var(--accent)'}">${fmt(totalSec)}</div>
      </div>
    </div>
    ${useCached ? '<div style="font-size:11px;color:var(--muted)">BayProX skipped (cached proxy). Fast run.</div>' : ''}
    ${tip}`;
}

function storeCol(key, role, value) {
  if (!colMap[key]) colMap[key] = {};
  colMap[key][role] = value;
}

// ── Config summary ───────────────────────────────────────────────────────────
function buildSummary() {
  // Dynamic TE rows from configured teRowData
  const teRows = (teRowData.length > 0 ? teRowData : [{col:'—',unit:''}])
    .map((r, i) => {
      const p = `te${i+1}`;
      const elem   = document.getElementById(`${p}_elem`)?.value || '—';
      const kdmn   = document.getElementById(`${p}_Kd_mn`)?.value || '—';
      const kdsd   = document.getElementById(`${p}_Kd_sd`)?.value || '—';
      const kp     = document.getElementById(`${p}_Kp`)?.value || '—';
      const unit   = r.unit || 'ppb';
      const col    = r.col  || '—';
      return [`TE${i+1} (${elem})`,
              `col: ${col} [${unit}] | Kd μ=${kdmn} σ=${kdsd} | Kp=${kp}`];
    });
  const modeLabel = analysisMode === 'full' ? '📐 Full quantification' : '📊 Semi-quantitative';
  const rows = [
    ['Station',          v('station_name')],
    ['Analysis mode',    modeLabel],
    ['Cal. age range',   v('calage_min') + ' – ' + v('calage_max') + ' yrs BP'],
    ...teRows,
    ['Cave temp',        v('temp_C') + ' °C'],
    ['Drip rate',        v('global_drip_rate') + ' drips/min'],
    ['Ca conc',          v('ca_conc') + ' ' + (document.getElementById('ca_unit')?.value||'ppb')],
    ['Realisations',     document.getElementById('generate_realisations')?.checked
                         ? v('n_realisations') : 'skipped'],
    ['Use cached proxy', document.getElementById('use_cached_proxy').checked ? 'Yes' : 'No'],
    ['V grid',           `V_max=${v('v_max')} drips/min, V_res=${v('v_res')}`],
  ];
  document.getElementById('summaryContent').innerHTML =
    rows.map(([k,val]) =>
      `<div style="display:flex;gap:16px;border-bottom:1px solid var(--border);padding:5px 0">
         <span style="color:var(--muted);min-width:160px">${k}</span>
         <span style="color:var(--texthi)">${val}</span>
       </div>`
    ).join('');
}

function v(id) { return document.getElementById(id)?.value ?? ''; }

// ── Run ──────────────────────────────────────────────────────────────────────
// ── Shared constants ────────────────────────────────────────────────────
const UNIT_FACTOR = {ppb:1, ppm:1000, 'ug/g':1000, 'mg/kg':1000};
// Ca fraction in CaCO₃ by mass: 40.078/100.087 ≈ 40% = 400,000 ppm
const CA_CALCITE_PPM = 400000;

// ── Pre-run sanity checks ────────────────────────────────────────────────

function runSanityChecks() {
  // Each warning: {level: 'error'|'warn', html: '...'}
  const warnings = [];
  const mode = analysisMode;

  const WARN_STYLE = {
    error: { bg: 'rgba(255,80,80,0.08)',   border: '#ff5050', color: '#ff5050' },
    warn:  { bg: 'rgba(255,160,50,0.08)',  border: '#ffa032', color: '#ffa032' },
  };
  function warnHtml(level, title, body) {
    const s = WARN_STYLE[level];
    return `<div style="margin-bottom:10px;padding:8px 10px;background:${s.bg};
                 border-left:3px solid ${s.border};border-radius:0 4px 4px 0">
              <strong style="color:${s.color}">${title}</strong><br>${body}
            </div>`;
  }

  teRowData.forEach((r, i) => {
    const n = i + 1;
    const p = `te${n}`;

    // 0. Depth unit mismatch info
    const daUnit = document.getElementById('da-depth-unit-sel')?.value || 'cm';
    if (i === 0 && daUnit !== teDepthUnit) {
      warnings.push({level: 'warn', html: warnHtml('warn',
        'Depth unit conversion active',
        `Depth/age CSV is in <strong>${daUnit}</strong> but TE data is in <strong>${teDepthUnit}</strong>. ` +
        `TE depths will be automatically converted to ${daUnit} before processing.`)});
    }

    // 1. Predicted vs observed speleothem concentration (model-correct)
    if (mode === 'full') {
      const aqVal    = parseFloat(document.getElementById(`${p}_aq_conc`)?.value) || 0;
      const aqUnit   = document.getElementById(`${p}_aq_unit`)?.value || 'ppb';
      const aqPPB    = aqVal * (UNIT_FACTOR[aqUnit] || 1);
      const caVal    = parseFloat(document.getElementById('ca_conc')?.value) || 0;
      const caUnit   = document.getElementById('ca_unit')?.value || 'ppb';
      const caPPB    = caVal * (UNIT_FACTOR[caUnit] || 1);
      const kdMn     = parseFloat(document.getElementById(`${p}_Kd_mn`)?.value) || 0;
      const xf       = parseFloat(document.getElementById(`${p}_F`)?.value) || 0;
      const xl       = parseFloat(document.getElementById(`${p}_labile`)?.value) || 0;
      const keVal    = parseFloat(document.getElementById(`${p}_K_e`)?.value) || 1;
      const kpVal    = parseFloat(document.getElementById(`${p}_Kp`)?.value);
      const kpEff    = (kpVal === -1 || !isFinite(kpVal))
                       ? (THEO_KP[document.getElementById(`${p}_elem`)?.value] || 1) : kpVal;
      const vRef     = getVRef();
      const tau      = 60.0 / vRef;
      const nS       = 1.0 - xf;
      const kd       = Math.exp(kdMn);
      const E1       = Math.exp(-kd * tau);
      const K0_ppm   = kpEff * (xf + xl) * (aqPPB / caPPB) * CA_CALCITE_PPM * keVal;
      const predPpm  = K0_ppm * (1 - nS * E1);
      const dataUnit = r.unit || 'ppb';
      const dataCol  = r.col;

      if (predPpm > 0 && teRawData[dataCol]) {
        const dataFactor = UNIT_FACTOR[dataUnit] || 1;
        const predNative = predPpm / (dataFactor / 1000);
        const rawVals    = (teRawData[dataCol] || []).map(Number).filter(isFinite);
        if (rawVals.length) {
          rawVals.sort((a,b)=>a-b);
          const obsNative = rawVals[Math.floor(rawVals.length/2)];
          const ratio     = obsNative / predNative;
          if (ratio > 5 || ratio < 0.2) {
            const dir  = ratio > 1 ? 'higher' : 'lower';
            warnings.push({ level: 'warn', html: warnHtml('warn',
              `TE${n} — prediction mismatch (${ratio.toFixed(1)}×) at ${vRef} drips/min`,
              `Kinetic prediction: <code>${predNative.toFixed(3)} ${dataUnit}</code>
               &nbsp;·&nbsp; Observed median: <code>${obsNative.toFixed(3)} ${dataUnit}</code>
               &nbsp;·&nbsp; ${ratio.toFixed(1)}× ${dir}.<br>
               <span style="color:var(--muted)">Adjust drip rate, Kd, or solution chemistry.</span>`)
            });
          }
        }
      }
    }

    // Note: aq_conc in ppb and data in ppm is normal (different physical units:
    // µg/L aqueous vs µg/g solid). No warning needed for this case.

    // 3. Missing aq_conc in full mode
    if (mode === 'full') {
      const aqVal = parseFloat(document.getElementById(`${p}_aq_conc`)?.value);
      if (!aqVal || aqVal <= 0) {
        warnings.push({ level: 'error', html: warnHtml('error',
          `TE${n} — aq_conc is zero or missing`,
          `Aqueous concentration must be > 0 in Full Quantification mode.`)
        });
      }
    }
  });

  // 4. Ca conc check
  if (mode === 'full') {
    const caVal  = parseFloat(document.getElementById('ca_conc')?.value) || 0;
    const caUnit = document.getElementById('ca_unit')?.value || 'ppb';
    const caMgL  = caVal * (UNIT_FACTOR[caUnit] || 1) / 1000;  // convert ppb → mg/L
    if (caMgL > 200) {
      warnings.push({ level: 'warn', html: warnHtml('warn',
        `Ca concentration unusually high (${caMgL.toFixed(0)} mg/L)`,
        `Typical cave drip water Ca is 20–100 mg/L. Check units and value.`)
      });
    }
  }

  const card = document.getElementById('sanity-card');
  const body = document.getElementById('sanity-body');
  if (!card || !body) return warnings.length === 0;
  if (warnings.length === 0) {
    card.style.display = 'none';
  } else {
    body.innerHTML = warnings.map(w => w.html).join('');
    card.style.display = '';
  }
  // true = ok to run (no error or warn-level warnings present)
  return warnings.filter(w => w.level === 'error' || w.level === 'warn').length === 0;
}

function startRun() {
  const ok = runSanityChecks();
  if (!ok) {
    const proceed = confirm(
      "Parameter warnings detected (see above).\n\n"
      + "The run may produce empty or NaN results.\n\n"
      + "Proceed anyway?"
    );
    if (!proceed) return;
  }
  const params = collectParams();
  if (!params) return;

  document.getElementById('runBtn').disabled = true;
  document.getElementById('resultsBtn').style.display = 'none';
  document.getElementById('logBox').innerHTML = '';
  setStatus('running', 'Running');
  runStartTime = Date.now();

  fetch('/run', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(params)
  })
  .then(r => r.json())
  .then(d => {
    if (d.error) { alert(d.error); resetRunBtn(); return; }
    pollTimer = setInterval(pollStatus, 1500);
  });
}

function collectParams() {
  const cm = colMap;
  return {
    station_name:   v('station_name'),
    calage_min:     v('calage_min'),
    calage_max:     v('calage_max'),
    // depth/age columns
    col_depth:      cm.depth_age?.depth   || '',
    col_age:        cm.depth_age?.age     || '',
    col_age_err:    cm.depth_age?.age_err || '',
    da_depth_unit:  document.getElementById('da-depth-unit-sel')?.value || 'cm',
    te_depth_unit:  teDepthUnit,
    // TE list — one entry per configured trace element row
    te_list: teRowData.map(r => ({col_depth: teDepthCol, col_proxy: r.col, unit: r.unit})),
    // TE param cards (te1_, te2_, te3_, ... collected dynamically)
    ...(() => {
      const out = {};
      for (let i = 1; i <= teRowData.length; i++) {
        const p = `te${i}`;
        ['elem','mol_wt','Kp','Kd_mn','Kd_sd','K_e','cal_pct','F','InertF','labile','aq_conc','aq_unit'].forEach(f => {
          const el = document.getElementById(`${p}_${f}`);
          if (el) out[`${p}_${f}`] = el.value;
        });
      }
      return out;
    })(),
    // isotope columns
    iso_col_depth:  cm.isotope1?.depth || '',
    iso_col_proxy:  cm.isotope1?.proxy || '',
    // cave
    temp_C:   v('temp_C'),
    global_drip_rate: v('global_drip_rate'),
    ca_conc:      v('ca_conc'),
    ca_unit:      document.getElementById('ca_unit')?.value || 'ppb',
    // Concentration priors (stochastic mode)
    ...(() => {
      const out = {};
      const caPr = concentrationPriors['ca'];
      if (caPr) {
        out['ca_prior_mu_ln']    = caPr.mu_ln;
        out['ca_prior_sigma_ln'] = caPr.sigma_ln;
        out['ca_prior_source']   = 'manual';
      }
      teRowData.forEach((_, i) => {
        const p  = `te${i+1}`;
        const pr = concentrationPriors[p];
        if (pr) {
          out[`${p}_prior_mu_ln`]    = pr.mu_ln;
          out[`${p}_prior_sigma_ln`] = pr.sigma_ln;
          out[`${p}_prior_source`]   = 'manual';
        }
      });
      return out;
    })(),
    // output
    n_realisations:       v('n_realisations'),
    rng_seed:             v('rng_seed'),
    generate_realisations: document.getElementById('generate_realisations')?.checked ?? true,
    analysis_mode:        analysisMode,
    v_max:                parseFloat(document.getElementById('v_max')?.value) || 100,
    v_res:                parseInt(document.getElementById('v_res')?.value) || 5000,
    hiatus_zones:         hiatusZones.filter(z => z.from !== z.to),
    outlier_win_size:     parseInt(document.getElementById('te-winsize')?.value) || 50,
    semi_anchor:          v('semi_anchor'),
    semi_ref_min:         v('semi_ref_min'),
    semi_ref_max:         v('semi_ref_max'),
    use_cached_proxy:     document.getElementById('use_cached_proxy').checked,
  };
}

function pollStatus() {
  fetch('/status').then(r => r.json()).then(s => {
    // progress bar
    document.getElementById('progressBar').style.width = s.progress + '%';
    document.getElementById('progressPct').textContent  = s.progress + '%';
    document.getElementById('stageLabel').textContent   = s.stage;

    // ETA — show elapsed until we have enough signal, then switch to ETA
    if (s.running && runStartTime) {
      const elapsed = (Date.now() - runStartTime) / 1000;
      const etaEl = document.getElementById('etaLabel');
      if (s.progress >= 5) {
        const remaining = (elapsed / s.progress) * (100 - s.progress);
        etaEl.textContent = 'ETA ~' + fmtTime(remaining);
      } else if (elapsed > 2) {
        etaEl.textContent = 'Elapsed: ' + fmtTime(elapsed);
      } else {
        etaEl.textContent = 'Starting…';
      }
    } else if (!s.running && s.done) {
      const elapsed = runStartTime ? (Date.now() - runStartTime) / 1000 : 0;
      if (elapsed > 0)
        document.getElementById('etaLabel').textContent = 'Completed in ' + fmtTime(elapsed);
    }

    // log
    const lb = document.getElementById('logBox');
    lb.innerHTML = s.log.map(line => {
      const cls = line.startsWith('✓') ? 'log-ok' : line.startsWith('ERROR') ? 'log-err' : '';
      return `<div class="log-line ${cls}">${escHtml(line)}</div>`;
    }).join('');
    lb.scrollTop = lb.scrollHeight;

    if (s.done) {
      clearInterval(pollTimer);
      document.getElementById('runBtn').disabled = false;
      if (s.error) {
        setStatus('error', 'Error');
      } else {
        setStatus('done', 'Complete');
        document.getElementById('resultsBtn').style.display = 'inline-flex';
        document.getElementById('badgeResults').style.display = 'inline';
        refreshDownloads(s.outputs);
        // Load hiatus zones from chart_data for results overlays
        fetch('/chart_data').then(r=>r.ok?r.json():null).then(d => {
          if (d && d.hiatus_zones) window._chartDataHiatusZones = d.hiatus_zones;
        }).catch(()=>{});
      }
    }
  });
}

function resetRunBtn() {
  document.getElementById('runBtn').disabled = false;
  setStatus('idle', 'Idle');
}

// ── Status pill ──────────────────────────────────────────────────────────────
function setStatus(state, label) {
  const pill = document.getElementById('statusPill');
  pill.className = 'status-pill ' + state;
  const dot = pill.querySelector('.dot');
  dot.className = 'dot' + (state === 'running' ? ' pulse' : '');
  pill.innerHTML = `<span class="dot${state==='running'?' pulse':''}"></span> ${label}`;
  // Toggle drip animation on sidebar image
  const navImg = document.querySelector('.nav-image-wrap');
  if (navImg) {
    if (state === 'running') navImg.classList.add('running');
    else navImg.classList.remove('running');
  }
}

// ── Chart ────────────────────────────────────────────────────────────────────
let sfChart = null;

function switchResultTab(name) {
  ['ts','sf','pdf','age'].forEach(t => {
    document.getElementById(`rtab-${t}`).style.display = t===name ? '' : 'none';
    const btn = document.getElementById(`tab-${t}`);
    btn.style.borderBottomColor = t===name ? 'var(--accent)' : 'transparent';
    btn.style.color = t===name ? 'var(--texthi)' : 'var(--muted)';
  });
  if (name==='sf') loadSFChart();
  if (name==='pdf') loadPdfHeatmap();
  if (name==='age') loadAgeModel();
}

function loadPdfHeatmap() {
  Promise.all([
    fetch('/data/pdf_heatmap.json').then(r => { if (!r.ok) throw new Error(r.status); return r.json(); }),
    fetch('/chart_data').then(r => r.ok ? r.json() : null).catch(() => null),
  ]).then(([hm, cd]) => {
    window._pdfHeatmapData = hm;
    window._pdfOverlayData = cd;  // percentiles for overlay
    if (cd && cd.hiatus_zones) window._chartDataHiatusZones = cd.hiatus_zones;
    setTimeout(renderPdfHeatmap, 80);
  }).catch(e => {
    console.warn('PDF heatmap load:', e);
    const c = document.getElementById('pdfHeatmapCanvas');
    if (c) {
      const ctx = c.getContext('2d');
      c.width = c.parentElement.clientWidth || 400;
      c.height = c.parentElement.clientHeight || 300;
      ctx.fillStyle = '#8fa4b5'; ctx.font = '13px Arial'; ctx.textAlign = 'center';
      ctx.fillText('No PDF data — run the model first', c.width/2, c.height/2);
    }
  });
}

function loadAgeModel() {
  fetch('/data/age_model.json')
    .then(r => { if (!r.ok) throw new Error(r.status); return r.json(); })
    .then(d => {
      window._ageModelData = d;
      setTimeout(renderAgeModelChart, 80);
    })
    .catch(e => {
      console.warn('Age model load:', e);
      const c = document.getElementById('ageModelChart');
      if (c) {
        const ctx = c.getContext('2d');
        c.width = c.parentElement.clientWidth || 400;
        c.height = c.parentElement.clientHeight || 300;
        ctx.fillStyle = '#8fa4b5'; ctx.font = '13px Arial'; ctx.textAlign = 'center';
        ctx.fillText('No age model data — run the model first', c.width/2, c.height/2);
      }
    });
}

function loadSFChart() {
  fetch('/chart_data').then(r => r.json()).then(d => {
    if (!d.sf) return;
    if (d.mode === 'semi') {
      const ctx = document.getElementById('sfChart').getContext('2d');
      if (sfChart) sfChart.destroy();
      ctx.font = '13px Arial, sans-serif';
      ctx.fillStyle = '#6e7f8d';
      ctx.textAlign = 'center';
      const h = ctx.canvas.height;
      ctx.fillText('Smart & Friedrich classification requires Full Quantification mode.',
                   ctx.canvas.width/2, h/2 - 10);
      ctx.fillText('Switch mode in Model Parameters to enable this chart.',
                   ctx.canvas.width/2, h/2 + 14);
      return;
    }
    const sf = d.sf;
    const ctx = document.getElementById('sfChart').getContext('2d');
    if (sfChart) sfChart.destroy();

    // Classification zones: [xMin, xMax, yMin, yMax, colour, label]
    const CV_THRESH  = 0.5;
    const Q_THRESH   = 5.0;   // drips per minute
    const zones = [
      // [cvLo, cvHi, qLo, qHi, rgba, name]
      [0,   CV_THRESH, 0,   Q_THRESH,  'rgba(91,155,213,0.13)',  'Seepage / percolation'],
      [CV_THRESH, 3,   0,   Q_THRESH,  'rgba(112,173,71,0.13)',  'Fracture / conduit'],
      [0,   CV_THRESH, Q_THRESH, 60,   'rgba(237,125,49,0.13)',  'Buffered overflow'],
      [CV_THRESH, 3,   Q_THRESH, 60,   'rgba(168,85,247,0.13)',  'Flood / conduit overflow'],
    ];

    // Zone background as custom plugin
    const zonePlugin = {
      id: 'sfZones',
      beforeDatasetsDraw(chart) {
        const {ctx:c, chartArea:{left,right,top,bottom}, scales:{x,y}} = chart;
        c.save();
        zones.forEach(([x0,x1,y0,y1,col]) => {
          const px0=x.getPixelForValue(x0), px1=x.getPixelForValue(x1);
          const py0=y.getPixelForValue(y0), py1=y.getPixelForValue(y1);
          c.fillStyle=col;
          c.fillRect(
            Math.max(left,Math.min(px0,px1)),
            Math.max(top, Math.min(py0,py1)),
            Math.abs(px1-px0), Math.abs(py1-py0)
          );
        });
        // Division lines
        c.strokeStyle='rgba(255,255,255,0.1)'; c.lineWidth=1; c.setLineDash([4,3]);
        const xMid=x.getPixelForValue(CV_THRESH);
        c.beginPath(); c.moveTo(xMid,top); c.lineTo(xMid,bottom); c.stroke();
        const yMid=y.getPixelForValue(Q_THRESH);
        c.beginPath(); c.moveTo(left,yMid); c.lineTo(right,yMid); c.stroke();
        c.restore();
      }
    };

    // Determine zone label for the point
    const inZone = zones.find(([x0,x1,y0,y1]) =>
      sf.cv>=x0 && sf.cv<x1 && sf.mean>=y0 && sf.mean<y1);
    const sfLabel = document.getElementById('sf-label');
    if (sfLabel) sfLabel.textContent = inZone
      ? `Classification: ${inZone[5]}` : '';

    sfChart = new Chart(ctx, {
      type: 'scatter',
      data: {
        datasets: [{
          label: 'This reconstruction',
          data: [{x: sf.cv, y: sf.mean}],
          backgroundColor: 'rgba(76,201,160,0.9)',
          pointRadius: 10,
          pointHoverRadius: 12,
        }]
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        plugins: {
          legend: {display: false},
          tooltip: {
            callbacks: {
              label: i => [
                `CV: ${sf.cv.toFixed(3)}  (range: ${sf.cv_lo.toFixed(3)}–${sf.cv_hi.toFixed(3)})`,
                `Mean drip rate: ${sf.mean.toFixed(2)} min⁻¹  (IQR: ${sf.mean_lo.toFixed(2)}–${sf.mean_hi.toFixed(2)})`,
              ]
            }
          }
        },
        scales: {
          x: {
            min: 0, max: 3,
            title: {display:true, text:'Coefficient of Variation (CV = σ/μ)',
                    color:'#8fa4b5', font:{size:11}},
            ticks: {color:'#8fa4b5', font:{size:10}},
            grid:  {color:'rgba(255,255,255,0.05)'},
          },
          y: {
            min: 0, max: 60,
            title: {display:true, text:'Mean drip rate (drips min⁻¹)',
                    color:'#8fa4b5', font:{size:11}},
            ticks: {color:'#8fa4b5', font:{size:10}},
            grid:  {color:'rgba(255,255,255,0.05)'},
          }
        }
      },
      plugins: [zonePlugin],
    });

    // Draw IQR error bars via afterDraw
    // (Chart.js scatter doesn't natively support error bars, draw manually)
    sfChart.options.animation = {onComplete: () => drawSFErrorBars(sfChart, sf)};
    sfChart.update();
  });
}

function drawSFErrorBars(chart, sf) {
  const {ctx:c, scales:{x,y}} = chart;
  const px = x.getPixelForValue(sf.cv);
  const py = y.getPixelForValue(sf.mean);
  const px0= x.getPixelForValue(sf.cv_lo),  px1= x.getPixelForValue(sf.cv_hi);
  const py0= y.getPixelForValue(sf.mean_lo), py1= y.getPixelForValue(sf.mean_hi);
  c.save();
  c.strokeStyle='rgba(76,201,160,0.6)'; c.lineWidth=2; c.setLineDash([]);
  // CV error bar (horizontal)
  c.beginPath(); c.moveTo(px0,py); c.lineTo(px1,py); c.stroke();
  c.beginPath(); c.moveTo(px0,py-5); c.lineTo(px0,py+5); c.stroke();
  c.beginPath(); c.moveTo(px1,py-5); c.lineTo(px1,py+5); c.stroke();
  // Mean error bar (vertical)
  c.beginPath(); c.moveTo(px,py0); c.lineTo(px,py1); c.stroke();
  c.beginPath(); c.moveTo(px-5,py0); c.lineTo(px+5,py0); c.stroke();
  c.beginPath(); c.moveTo(px-5,py1); c.lineTo(px+5,py1); c.stroke();
  c.restore();
}

// ── Y-axis mode: 'drip' (drips/min) or 'tau' (seconds) ──────────────────
let yAxisMode = 'drip';

function setYMode(mode) {
  yAxisMode = mode;
  document.getElementById('ymode-drip').classList.toggle('active', mode === 'drip');
  document.getElementById('ymode-tau').classList.toggle('active', mode === 'tau');
  const title = document.getElementById('ts-card-title');
  if (title) title.textContent = mode === 'drip' ? '📈 Drip Rate vs Time' : '📈 Residence Time (τ) vs Time';
  // Re-render both charts
  loadChart();
  if (window._pdfHeatmapData) renderPdfHeatmap();
}

// Convert drip rate (drips/min) to τ (seconds): τ = 60/d
function dripToTau(v) {
  if (v === null || v === undefined || !isFinite(v) || v <= 0) return null;
  return 60.0 / v;
}

function loadChart() {
  fetch('/chart_data').then(r => r.json()).then(d => {
    if (d.error) return;
    const ctx = document.getElementById('dripChart').getContext('2d');
    if (chart) chart.destroy();

    // Build {x, y} point arrays, using null for NaN to create gaps
    // When yAxisMode='tau', convert drip rate → τ = 60/d
    const useTau = yAxisMode === 'tau';
    const mkPts = arr => d.age.map((a, i) => {
      const v = arr[i];
      if (v === null || v === undefined || !isFinite(v)) return {x: a, y: null};
      return {x: a, y: useTau ? dripToTau(v) : v};
    });

    // Hiatus zone overlay plugin
    const tsHiatusPlugin = {
      id: 'tsHiatus',
      afterDraw(chart) {
        const zones = d.hiatus_zones || [];
        if (!zones.length) return;
        const {ctx: c, chartArea: {left,right,top,bottom}, scales: {x}} = chart;
        c.save();
        zones.forEach(z => {
          const xF = x.getPixelForValue(z.from), xT = x.getPixelForValue(z.to);
          const xL = Math.min(xF, xT), xR = Math.max(xF, xT);
          c.fillStyle = 'rgba(255,80,80,0.1)';
          c.fillRect(xL, top, xR-xL, bottom-top);
          c.strokeStyle = 'rgba(255,80,80,0.3)';
          c.setLineDash([2,2]); c.lineWidth = 1;
          c.beginPath(); c.moveTo(xL,top); c.lineTo(xL,bottom); c.stroke();
          c.beginPath(); c.moveTo(xR,top); c.lineTo(xR,bottom); c.stroke();
          c.fillStyle = 'rgba(255,80,80,0.5)'; c.font = '10px Arial';
          c.textAlign = 'center';
          c.fillText('hiatus', (xL+xR)/2, top + 14);
        });
        c.restore();
      }
    };

    chart = new Chart(ctx, {
      type: 'line',
      data: {
        datasets: [
          { label: 'pc95\u2013pc05 range', data: mkPts(d.pc95), fill: '-1',
            backgroundColor: 'rgba(76,201,160,0.08)', borderColor: 'transparent',
            pointRadius: 0, tension: 0.3, spanGaps: false },
          { label: 'pc05', data: mkPts(d.pc05), borderColor: 'rgba(76,201,160,0.25)',
            borderWidth: 1, pointRadius: 0, fill: false, tension: 0.3, spanGaps: false },
          { label: 'IQR (25\u201375)', data: mkPts(d.pc75), fill: '+1',
            backgroundColor: 'rgba(76,201,160,0.15)', borderColor: 'transparent',
            pointRadius: 0, tension: 0.3, spanGaps: false },
          { label: 'pc25', data: mkPts(d.pc25), borderColor: 'rgba(76,201,160,0.35)',
            borderWidth: 1, pointRadius: 0, fill: false, tension: 0.3, spanGaps: false },
          { label: 'Median', data: mkPts(d.pc50), borderColor: '#4cc9a0',
            borderWidth: 2, pointRadius: 0, fill: false, tension: 0.3, spanGaps: false },
          ...(d.kd_lo && d.kd_hi ? [
            { label: 'Kd \u22121\u03c3', data: mkPts(d.kd_lo),
              borderColor: 'rgba(200,180,80,0.45)', borderWidth: 1, borderDash: [3,3],
              pointRadius: 0, fill: false, tension: 0.3, spanGaps: false },
            { label: 'Kd +1\u03c3', data: mkPts(d.kd_hi),
              borderColor: 'rgba(200,180,80,0.45)', borderWidth: 1, borderDash: [3,3],
              pointRadius: 0, fill: false, tension: 0.3, spanGaps: false },
          ] : []),
        ]
      },
      options: {
        responsive: true, maintainAspectRatio: false, animation: false,
        interaction: { mode: 'index', intersect: false },
        plugins: {
          legend: { labels: { color: '#8fa4b5', font: { size: 11 } } },
          decimation: { enabled: true, algorithm: 'lttb', samples: 800 },
        },
        scales: {
          x: { type: 'linear', reverse: true,
               ticks: { color: '#8fa4b5', font: { size: 10 } },
               grid: { color: 'rgba(42,52,65,0.6)' },
               title: { display: true, text: 'Age (yrs BP)', color: '#8fa4b5', font: { size: 11 } } },
          y: { type: 'linear',
               ticks: { color: '#8fa4b5', font: { size: 10 } },
               grid: { color: 'rgba(42,52,65,0.6)' },
               title: { display: true,
                 text: d.mode === 'semi' ? 'Relative drip rate (% of reference)'
                       : useTau ? 'Residence time \u03c4 (s)' : 'Drip rate (min\u207b\u00b9)',
                 color: '#8fa4b5', font: { size: 11 } } }
        }
      },
      plugins: [tsHiatusPlugin],
    });
  }).catch(e => console.error('loadChart error:', e));
}

// ── Downloads ────────────────────────────────────────────────────────────────
// ── Colourmap helpers ─────────────────────────────────────────────────────
const CMAPS = {
  viridis:  [[68,1,84],[72,36,117],[65,68,135],[53,95,141],[42,120,142],[33,145,140],[34,168,132],[53,191,111],[94,211,67],[163,222,21],[253,231,37]],
  inferno:  [[0,0,4],[22,11,57],[66,10,104],[106,23,110],[147,38,103],[188,55,84],[221,81,58],[243,118,27],[252,166,10],[246,215,70],[252,255,164]],
  plasma:   [[13,8,135],[75,3,161],[126,3,168],[168,34,150],[198,62,118],[221,95,84],[236,130,55],[245,167,27],[246,207,18],[231,244,40],[240,249,33]],
  magma:    [[0,0,4],[18,13,49],[51,16,104],[89,26,120],[128,40,118],[165,63,111],[199,90,103],[227,126,101],[248,170,109],[254,216,144],[252,253,191]],
  cividis:  [[0,34,78],[0,57,100],[38,77,108],[73,96,111],[107,114,117],[140,132,120],[172,151,112],[202,172,93],[228,195,64],[248,222,30],[255,255,0]],
  turbo:    [[48,18,59],[69,91,205],[32,158,245],[18,209,192],[66,243,114],[147,254,57],[220,237,30],[255,193,22],[250,131,13],[217,65,8],[122,4,3]],
  greens:   [[247,252,245],[229,245,224],[199,233,192],[161,217,155],[116,196,118],[65,171,93],[35,139,69],[0,109,44],[0,68,27],[0,48,18],[0,32,10]],
  blues:    [[247,251,255],[222,235,247],[198,219,239],[158,202,225],[107,174,214],[66,146,198],[33,113,181],[8,81,156],[8,48,107],[3,32,76],[1,18,50]],
  reds:     [[255,245,240],[254,224,210],[252,187,161],[252,146,114],[251,106,74],[239,59,44],[203,24,29],[165,15,21],[103,0,13],[70,0,8],[45,0,4]],
};

function sampleCmap(name, t) {
  const c = CMAPS[name] || CMAPS.viridis;
  const n = c.length - 1;
  const i = Math.min(Math.floor(t * n), n - 1);
  const f = t * n - i;
  return c[i].map((v, k) => Math.round(v + f * (c[i+1][k] - v)));
}

// ── PDF Heatmap rendering ────────────────────────────────────────────────
let _pdfRetries = 0;
function renderPdfHeatmap() {
  const canvas = document.getElementById('pdfHeatmapCanvas');
  if (!canvas || !window._pdfHeatmapData) return;
  const parentW = canvas.parentElement.clientWidth;
  const parentH = canvas.parentElement.clientHeight;
  if (!parentW || !parentH) {
    if (_pdfRetries++ < 10) setTimeout(renderPdfHeatmap, 100);
    return;
  }
  _pdfRetries = 0;

  const {V_pdf, V_span, ages} = window._pdfHeatmapData;
  const nv = V_pdf.length, nt = V_pdf[0].length;
  const cmap = document.getElementById('pdf-cmap')?.value || 'viridis';
  const useLog = document.getElementById('pdf-log')?.checked || false;

  // Normalize each timestep column to [0,1] by its own max (per-column)
  // This ensures every timestep shows its PDF shape at full contrast
  const normPdf = [];
  for (let v = 0; v < nv; v++) normPdf.push(new Array(nt).fill(0));
  for (let t = 0; t < nt; t++) {
    let colMax = 0;
    for (let v = 0; v < nv; v++) {
      const val = V_pdf[v][t] || 0;
      if (val > colMax) colMax = val;
    }
    if (colMax > 0) {
      for (let v = 0; v < nv; v++) normPdf[v][t] = (V_pdf[v][t] || 0) / colMax;
    }
  }

  let maxD = 1.0;  // already normalised to [0,1] per column
  const data = [];
  for (let v = 0; v < nv; v++) {
    const row = [];
    for (let t = 0; t < nt; t++) {
      let val = normPdf[v][t];
      if (useLog && val > 0) val = Math.log10(val * 9 + 1);  // log scale mapped to [0,1]
      row.push(val);
    }
    data.push(row);
  }

  const dpr = window.devicePixelRatio || 1;
  const W = parentW;
  const H = parentH;
  canvas.width = W * dpr; canvas.height = H * dpr;
  canvas.style.width = W + 'px'; canvas.style.height = H + 'px';
  const ctx = canvas.getContext('2d');
  ctx.scale(dpr, dpr);

  const ML = 60, MR = 40, MT = 10, MB = 35;
  const pw = W - ML - MR, ph = H - MT - MB;

  const tmpCanvas = document.createElement('canvas');
  tmpCanvas.width = nt; tmpCanvas.height = nv;
  const tmpCtx = tmpCanvas.getContext('2d');
  const img = tmpCtx.createImageData(nt, nv);
  for (let v = 0; v < nv; v++) {
    for (let t = 0; t < nt; t++) {
      // Reverse x: oldest on left, youngest on right (matches time series)
      const val = data[nv - 1 - v][nt - 1 - t];
      const norm = isFinite(val) && maxD > 0 ? Math.max(0, Math.min(1, val / maxD)) : 0;
      const [r, g, b] = sampleCmap(cmap, norm);
      const idx = (v * nt + t) * 4;
      img.data[idx] = r; img.data[idx+1] = g; img.data[idx+2] = b; img.data[idx+3] = 255;
    }
  }
  tmpCtx.putImageData(img, 0, 0);
  ctx.imageSmoothingEnabled = false;
  ctx.drawImage(tmpCanvas, ML, MT, pw, ph);

  // Axes (reversed: oldest left, youngest right)
  ctx.strokeStyle = '#8fa4b5'; ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(ML, MT); ctx.lineTo(ML, MT+ph); ctx.lineTo(ML+pw, MT+ph); ctx.stroke();
  ctx.fillStyle = '#8fa4b5'; ctx.font = '11px Arial, sans-serif'; ctx.textAlign = 'center';
  const ageMin = ages[0], ageMax = ages[nt-1];
  for (let i = 0; i <= 6; i++) {
    // frac=0 → left → oldest; frac=1 → right → youngest
    const frac = i / 6, age = ageMax - frac * (ageMax - ageMin), px = ML + frac * pw;
    ctx.fillText(Math.round(age).toString(), px, MT + ph + 15);
    ctx.beginPath(); ctx.moveTo(px, MT+ph); ctx.lineTo(px, MT+ph+4); ctx.stroke();
  }
  ctx.fillText('Age (yrs BP)', ML + pw/2, MT + ph + 30);
  const vMin = V_span[0], vMax = V_span[nv-1];
  const useTau = yAxisMode === 'tau';
  ctx.textAlign = 'right';
  for (let i = 0; i <= 5; i++) {
    const frac = i / 5, v = vMin + frac * (vMax - vMin), py = MT + ph - frac * ph;
    const label = useTau ? (v > 0 ? (60/v).toFixed(1) : '\u221e') : v.toFixed(1);
    ctx.fillText(label, ML - 5, py + 4);
    ctx.beginPath(); ctx.moveTo(ML, py); ctx.lineTo(ML-3, py); ctx.stroke();
  }
  ctx.save(); ctx.translate(14, MT + ph/2); ctx.rotate(-Math.PI/2);
  ctx.textAlign = 'center';
  ctx.fillText(useTau ? 'Residence time \u03c4 (s)' : 'Drip rate (min\u207b\u00b9)', 0, 0);
  ctx.restore();

  // ── Percentile overlay ────────────────────────────────────────────
  const showOverlay = document.getElementById('pdf-overlay')?.checked;
  const cd = window._pdfOverlayData;
  if (showOverlay && cd && cd.age && cd.pc50) {
    // Map age → pixel x (reversed: high age = left)
    const ageToPx = a => ML + (1 - (a - ageMin) / (ageMax - ageMin)) * pw;
    // Map drip rate → pixel y
    const vToPy = v => MT + ph - ((v - vMin) / (vMax - vMin)) * ph;

    // Subsample to ~500 points for smooth lines
    const step = Math.max(1, Math.floor(cd.age.length / 500));

    function drawLine(arr, color, width, dash) {
      ctx.beginPath(); ctx.strokeStyle = color; ctx.lineWidth = width;
      ctx.setLineDash(dash || []);
      let started = false;
      for (let i = 0; i < cd.age.length; i += step) {
        const v = arr[i];
        if (v === null || v === undefined || !isFinite(v)) { started = false; continue; }
        const px = ageToPx(cd.age[i]), py = vToPy(v);
        // Clip to plot area
        if (py < MT || py > MT + ph) { started = false; continue; }
        if (!started) { ctx.moveTo(px, py); started = true; }
        else ctx.lineTo(px, py);
      }
      ctx.stroke();
      ctx.setLineDash([]);
    }

    // IQR fill
    ctx.beginPath();
    ctx.fillStyle = 'rgba(255,255,255,0.08)';
    let pts25 = [], pts75 = [];
    for (let i = 0; i < cd.age.length; i += step) {
      const v25 = cd.pc25[i], v75 = cd.pc75[i];
      if (v25 != null && v75 != null && isFinite(v25) && isFinite(v75)) {
        pts25.push({px: ageToPx(cd.age[i]), py: vToPy(v25)});
        pts75.push({px: ageToPx(cd.age[i]), py: vToPy(v75)});
      }
    }
    if (pts25.length > 2) {
      ctx.beginPath();
      pts25.forEach((p, i) => i === 0 ? ctx.moveTo(p.px, p.py) : ctx.lineTo(p.px, p.py));
      [...pts75].reverse().forEach(p => ctx.lineTo(p.px, p.py));
      ctx.closePath();
      ctx.fill();
    }

    // Lines: pc05, pc25, median, pc75, pc95
    if (cd.pc05) drawLine(cd.pc05, 'rgba(255,255,255,0.2)', 0.8, [2,2]);
    if (cd.pc25) drawLine(cd.pc25, 'rgba(255,255,255,0.3)', 0.8);
    if (cd.pc50) drawLine(cd.pc50, 'rgba(255,255,255,0.9)', 1.5);
    if (cd.pc75) drawLine(cd.pc75, 'rgba(255,255,255,0.3)', 0.8);
    if (cd.pc95) drawLine(cd.pc95, 'rgba(255,255,255,0.2)', 0.8, [2,2]);

    // Legend
    ctx.fillStyle = 'rgba(255,255,255,0.7)'; ctx.font = '9px Arial'; ctx.textAlign = 'left';
    const lx = ML + 8, ly = MT + 12;
    ctx.fillStyle = 'rgba(255,255,255,0.9)';
    ctx.fillRect(lx, ly, 2, 8); ctx.fillText('Median', lx + 6, ly + 8);
    ctx.fillStyle = 'rgba(255,255,255,0.3)';
    ctx.fillRect(lx, ly + 12, 2, 8); ctx.fillText('IQR', lx + 6, ly + 20);
    ctx.fillStyle = 'rgba(255,255,255,0.2)';
    ctx.fillRect(lx, ly + 24, 2, 8); ctx.fillText('5–95%', lx + 6, ly + 32);
  }

  // Colorbar (right side)
  const cbW = 12, cbH = ph, cbX = ML + pw + 4, cbY = MT;
  for (let j = 0; j < cbH; j++) {
    const norm = 1 - j / cbH;  // top=1, bottom=0
    const [r, g, b] = sampleCmap(cmap, norm);
    ctx.fillStyle = `rgb(${r},${g},${b})`;
    ctx.fillRect(cbX, cbY + j, cbW, 1);
  }
  ctx.strokeStyle = '#8fa4b5'; ctx.lineWidth = 0.5;
  ctx.strokeRect(cbX, cbY, cbW, cbH);
  ctx.fillStyle = '#8fa4b5'; ctx.font = '9px Arial, sans-serif'; ctx.textAlign = 'left';
  ctx.fillText('1.0', cbX + cbW + 3, cbY + 9);
  ctx.fillText('0.5', cbX + cbW + 3, cbY + cbH/2 + 3);
  ctx.fillText('0.0', cbX + cbW + 3, cbY + cbH - 1);
  ctx.fillText('P', cbX + cbW + 3, cbY - 4);

  const xr = document.getElementById('pdf-xrange');
  const yr = document.getElementById('pdf-yrange');
  if (xr) xr.textContent = `${nt} timesteps \u00b7 column-normalised PDF (0\u20131)`;
  if (yr) yr.textContent = `${nv} V bins \u00b7 ${cmap}${useLog ? ' (log)' : ''}`;
}

// ── Age model chart ──────────────────────────────────────────────────────
let ageModelChart = null;

let _ageRetries = 0;
function renderAgeModelChart() {
  const canvas = document.getElementById('ageModelChart');
  if (!canvas || !window._ageModelData) return;
  if (!canvas.parentElement.clientWidth) {
    if (_ageRetries++ < 10) setTimeout(renderAgeModelChart, 100);
    return;
  }
  _ageRetries = 0;
  const {depth, age_median, age_lo, age_hi, dated_depth, dated_age, dated_err} = window._ageModelData;
  if (ageModelChart) ageModelChart.destroy();

  // Subsample to ~500 points for Chart.js performance
  const N = depth.length;
  const step = Math.max(1, Math.floor(N / 500));
  const subD = [], subA = [];
  for (let i = 0; i < N; i += step) {
    subD.push(depth[i]);
    subA.push(age_median[i]);
  }
  // Always include last point
  if (subD[subD.length-1] !== depth[N-1]) {
    subD.push(depth[N-1]);
    subA.push(age_median[N-1]);
  }

  const datasets = [];
  // Median age model (subsampled)
  datasets.push({
    label: 'BayProX median',
    data: subD.map((d, i) => ({x: d, y: subA[i]})),
    type: 'line',
    borderColor: 'rgba(76,201,160,0.9)', borderWidth: 2,
    pointRadius: 0, fill: false,
  });
  // Dated points with error bars
  if (dated_depth && dated_depth.length) {
    const errBars = [];
    for (let i = 0; i < dated_depth.length; i++) {
      const e = dated_err ? dated_err[i] : 0;
      if (e > 0) errBars.push({x:dated_depth[i],y:dated_age[i]-e},{x:dated_depth[i],y:dated_age[i]+e},{x:NaN,y:NaN});
    }
    if (errBars.length) datasets.push({
      label: '±1σ error', data: errBars, type: 'line',
      borderColor: 'rgba(247,164,64,0.4)', borderWidth: 1.5,
      pointRadius: 0, fill: false, spanGaps: false,
    });
    datasets.push({
      label: 'Dated points',
      data: dated_depth.map((d,i) => ({x:d, y:dated_age[i]})),
      type: 'scatter',
      backgroundColor: 'rgba(247,164,64,0.9)', pointRadius: 5,
    });
  }

  // Hiatus zone overlay plugin for results age model
  const ageHiatusPlugin = {
    id: 'ageModelHiatus',
    afterDraw(chart) {
      // Try loading hiatus zones from chart_data if available
      const zones = window._chartDataHiatusZones || [];
      if (!zones.length) return;
      const {ctx: c, chartArea: {left,right,top,bottom}, scales: {y}} = chart;
      c.save();
      zones.forEach(z => {
        const yF = y.getPixelForValue(z.from), yT = y.getPixelForValue(z.to);
        const yTop = Math.min(yF, yT), yBot = Math.max(yF, yT);
        c.fillStyle = 'rgba(255,80,80,0.12)';
        c.fillRect(left, yTop, right-left, yBot-yTop);
        c.strokeStyle = 'rgba(255,80,80,0.4)';
        c.setLineDash([2,2]); c.lineWidth = 1;
        c.beginPath(); c.moveTo(left,yTop); c.lineTo(right,yTop); c.stroke();
        c.beginPath(); c.moveTo(left,yBot); c.lineTo(right,yBot); c.stroke();
        c.fillStyle = 'rgba(255,80,80,0.6)'; c.font = '9px Arial, sans-serif';
        c.textAlign = 'center';
        c.fillText('hiatus', (left+right)/2, (yTop+yBot)/2 + 4);
      });
      c.restore();
    }
  };

  ageModelChart = new Chart(canvas, {
    data: {datasets},
    options: {
      animation: false, responsive: true, maintainAspectRatio: false,
      plugins: { legend: { labels: { color:'#cdd9e5', font:{size:10}, filter: i => !i.text.startsWith('_') } } },
      scales: {
        x: { type: 'linear', title:{display:true, text:'Depth (' + teDepthUnit + ')', color:'#8fa4b5', font:{size:11}},
             grid:{color:'rgba(255,255,255,0.05)'}, ticks:{color:'#8fa4b5',font:{size:10}} },
        y: { type: 'linear', title:{display:true, text:'Age (yrs BP)', color:'#8fa4b5', font:{size:11}},
             grid:{color:'rgba(255,255,255,0.05)'}, ticks:{color:'#8fa4b5',font:{size:10}} },
      }
    },
    plugins: [ageHiatusPlugin],
  });
}

function refreshDownloads(files) {
  if (!files) {
    fetch('/status').then(r=>r.json()).then(s => refreshDownloads(s.outputs));
    return;
  }
  const list = document.getElementById('downloadList');
  if (!files || files.length === 0) {
    list.innerHTML = '<div style="color:var(--muted);font-size:12px">No outputs yet — run the model first.</div>';
    return;
  }
  const descs = {
    'drip_rate_summary.csv':      'Percentile summary (pc05–pc95) at each time step',
    'drip_rate_realisations.csv': 'Full ensemble of realisations for RQA analysis',
    'age_model.json':             'BayProX age-depth model (browser chart data)',
    'age_model.csv':              'Age-depth model as CSV (depth, age_yBP)',
    'pdf_heatmap.json':           'Drip rate PDF matrix (subsampled for browser)',
    'input_summary.csv':          'All input parameters for this run (verification)',
    'run_id.txt':                 'Unique run identifier (links outputs to inputs)',
    'ProxyRecord.pkl':            'BayProX proxy record (reusable cache)',
  };
  list.innerHTML = files.map(fn => `
    <div class="output-item">
      <span class="oi-icon">📄</span>
      <div>
        <div class="oi-name">${fn}</div>
        <div class="oi-desc">${descs[fn] || ''}</div>
      </div>
      <a href="/download/${fn}" class="btn btn-ghost" style="padding:6px 14px;font-size:11px">Download</a>
    </div>`).join('');
}


// ── Theoretical Kp display (Wang & Xu 2001) ────────────────────────────────────
// When Kp = -1, compute and show theoretical value in adjacent read-only field.
// Note: for Co and Ni the model uses empirical (Lindeman) values; for other elements
// the Wang & Xu lattice-strain formula is used.
const THEO_KP = { 'Co': 4.4, 'Ni': 1.1, 'Cu': 44 };  // Lindeman et al. GCA 2022

// ── Reactive parameter hints ────────────────────────────────────────────

// Returns the raw median of proxy data in its native units (no conversion).
// Used for model comparison — the forward model works in data-native units.
// The unit selector is informational only; the data is NOT converted before
// entering the PDist / BayProX pipeline.
function getCalPct(prefix) {
  // Read calibration window percentage from slider (default 100 = full record)
  const el = document.getElementById(prefix + '_cal_pct');
  return el ? parseInt(el.value) || 100 : 100;
}

function getObsMedianNative(teIdx, calPct) {
  // Returns the median of the proxy column, optionally restricted to the
  // shallowest (most recent) calPct% of the depth range.
  const row = teRowData[teIdx];
  if (!row || !teRawData[row.col]) return null;
  const proxyVals = (teRawData[row.col]||[]).map(v=>parseFloat(v));
  const depthVals = teRawDepth;

  if (!proxyVals.length) return null;
  calPct = calPct || 100;

  if (calPct >= 100 || !depthVals || depthVals.length !== proxyVals.length) {
    // Full record — original behaviour
    const vals = proxyVals.filter(isFinite);
    if (!vals.length) return null;
    vals.sort((a,b)=>a-b);
    return vals[Math.floor(vals.length/2)];
  }

  // Depth-windowed: shallowest calPct% of the depth range
  const finiteDepths = depthVals.filter(isFinite);
  if (!finiteDepths.length) return null;
  const depMin = finiteDepths.reduce((a,b) => a < b ? a : b, Infinity);
  const depMax = finiteDepths.reduce((a,b) => a > b ? a : b, -Infinity);
  const depRange = depMax - depMin;
  if (depRange <= 0) return null;

  // "Most recent" = shallowest = smallest depth values
  const depThreshold = depMin + depRange * (calPct / 100);
  const windowVals = [];
  for (let i = 0; i < proxyVals.length; i++) {
    if (isFinite(proxyVals[i]) && isFinite(depthVals[i]) && depthVals[i] <= depThreshold) {
      windowVals.push(proxyVals[i]);
    }
  }
  if (!windowVals.length) return null;
  windowVals.sort((a,b)=>a-b);
  return windowVals[Math.floor(windowVals.length/2)];
}

function fmtPPB(v) {
  if (!isFinite(v)||v<=0) return null;
  if (v>=1000) return (v/1000).toFixed(3)+' ppm';
  if (v>=1)    return v.toFixed(3)+' ppb';
  return (v*1000).toFixed(1)+' ppt';
}

// ── Concentration prior state and fitting ────────────────────────────────
const concentrationPriors = {};   // keyed by 'ca', 'te1', 'te2', etc.

function fitLognormal(values_ppb) {
  const v = values_ppb.filter(x => isFinite(x) && x > 0);
  if (v.length < 3) return null;
  const lnv = v.map(Math.log);
  const mu = lnv.reduce((a,b)=>a+b,0) / lnv.length;
  const sigma = Math.sqrt(lnv.map(x=>(x-mu)**2).reduce((a,b)=>a+b,0) / (lnv.length-1));
  return { mu_ln: mu, sigma_ln: sigma };
}

function fitLognormalFromMeanSd(mean_ppb, sd_ppb) {
  if (mean_ppb <= 0 || sd_ppb <= 0) return null;
  sd_ppb = Math.max(sd_ppb, mean_ppb * 0.01);
  const sigma2 = Math.log(1 + (sd_ppb/mean_ppb)**2);
  const mu = Math.log(mean_ppb) - 0.5*sigma2;
  return { mu_ln: mu, sigma_ln: Math.sqrt(sigma2) };
}

function priorSummaryHTML(prior, label) {
  if (!prior) return '';
  const mean_ppb = Math.exp(prior.mu_ln + 0.5*prior.sigma_ln**2);
  const cv = Math.sqrt(Math.exp(prior.sigma_ln**2) - 1) * 100;
  return `<strong>${label} prior:</strong>
    µ_ln=${prior.mu_ln.toFixed(3)}, σ_ln=${prior.sigma_ln.toFixed(3)} —
    mean=${mean_ppb.toFixed(2)} ppb, CV=${cv.toFixed(1)}% (lognormal)`;
}

function describeDistribution(values_ppb, mu_ln, sigma_ln) {
  const cv = Math.sqrt(Math.exp(sigma_ln**2) - 1);
  let shape, fit, warning = null;
  if (sigma_ln < 0.15) {
    shape = 'approximately symmetric — low variability';
    fit = 'Lognormal fitted but near-normal behaviour.';
  } else if (sigma_ln < 0.4) {
    shape = 'mildly right-skewed';
    fit = 'Lognormal is a good fit. Moderate variability.';
  } else if (sigma_ln < 0.8) {
    shape = 'moderately right-skewed — typical for trace elements';
    fit = 'Lognormal is appropriate. Long right tail reflects episodic higher-concentration events.';
  } else if (sigma_ln < 1.3) {
    shape = 'strongly right-skewed — high variability';
    fit = 'High sigma_ln will substantially widen the posterior. Check for outliers.';
    warning = 'CV > 80%. Review raw data for outliers.';
  } else {
    shape = 'very strongly skewed — possible outliers';
    fit = 'Extreme skew — posterior reliability may be reduced.';
    warning = '\u26a0 sigma_ln > 1.3. Consider removing outliers or using manual mean + SD.';
  }
  return { shape, fit, warning, cv };
}

function renderConcentrationPlot(values_ppb, prior, canvasId, descId, label) {
  const v = values_ppb.filter(x => isFinite(x) && x > 0).sort((a,b) => a-b);
  if (v.length < 3) return;
  const { mu_ln, sigma_ln } = prior;
  const { shape, fit, warning, cv } = describeDistribution(v, mu_ln, sigma_ln);
  // Log-spaced histogram bins
  const N_BINS = Math.min(20, Math.ceil(Math.sqrt(v.length)));
  const logMin = Math.log(v[0]*0.95), logMax = Math.log(v[v.length-1]*1.05);
  const binEdges = Array.from({length:N_BINS+1},(_, i)=>Math.exp(logMin+i*(logMax-logMin)/N_BINS));
  const counts = new Array(N_BINS).fill(0);
  v.forEach(x => { const i = binEdges.findIndex((e,j)=>j<N_BINS&&x>=binEdges[j]&&x<binEdges[j+1]);
    if (i>=0) counts[i]++; });
  const binWidths = binEdges.slice(1).map((e,i)=>e-binEdges[i]);
  const densities = counts.map((c,i)=>c/(v.length*binWidths[i]));
  const binCentres = binEdges.slice(1).map((e,i)=>(e+binEdges[i])/2);
  // Lognormal PDF curve
  const N_CURVE = 120;
  const xCurve = Array.from({length:N_CURVE},(_, i)=>Math.exp(logMin+i*(logMax-logMin)/(N_CURVE-1)));
  const pdfCurve = xCurve.map(x => {
    const lx = Math.log(x);
    return Math.exp(-0.5*((lx-mu_ln)/sigma_ln)**2) / (x*sigma_ln*Math.sqrt(2*Math.PI));
  });
  const canvas = document.getElementById(canvasId);
  if (!canvas) return;
  const existing = Chart.getChart(canvas);
  if (existing) existing.destroy();
  const chart = new Chart(canvas, {
    data: { datasets: [
      { type:'bar', label:'Observed',
        data: binCentres.map((x,i)=>({x, y:densities[i]})),
        backgroundColor:'rgba(76,201,160,0.35)',
        borderColor:'rgba(76,201,160,0.8)', borderWidth:1,
        barPercentage:0.95, categoryPercentage:1.0 },
      { type:'line', label:'Lognormal fit',
        data: xCurve.map((x,i)=>({x, y:pdfCurve[i]})),
        borderColor:'#f7a440', borderWidth:2,
        pointRadius:0, tension:0.4, fill:false },
    ]},
    options: {
      scales: {
        x: { type:'linear',
             title:{display:true, text:label+' (ppb)', color:'#8fa4b5', font:{size:10}},
             ticks:{color:'#8fa4b5', font:{size:9},
               callback: v=>v>=1000?(v/1000).toFixed(1)+'k':v.toFixed(1)} },
        y: { title:{display:true, text:'Density', color:'#8fa4b5', font:{size:10}},
             ticks:{color:'#8fa4b5', font:{size:9}}, beginAtZero:true },
      },
      plugins: { legend:{labels:{color:'#cdd9e5',font:{size:10},boxWidth:12}} },
      animation:false, responsive:true, maintainAspectRatio:false,
    }
  });
  window[canvasId+'_chart'] = chart;
  window[canvasId+'_logScale'] = false;
  // Description block
  const desc = document.getElementById(descId);
  if (!desc) return;
  const mean_ppb = Math.exp(mu_ln + 0.5*sigma_ln**2);
  const med_ppb = Math.exp(mu_ln);
  desc.innerHTML = `
    <div style="display:grid;grid-template-columns:repeat(4,1fr);gap:6px;margin-bottom:8px">
      <div><span style="color:var(--muted)">n</span><br><strong>${v.length}</strong></div>
      <div><span style="color:var(--muted)">median</span><br><strong>${med_ppb.toFixed(2)} ppb</strong></div>
      <div><span style="color:var(--muted)">mean</span><br><strong>${mean_ppb.toFixed(2)} ppb</strong></div>
      <div><span style="color:var(--muted)">CV</span><br><strong>${(cv*100).toFixed(0)}%</strong></div>
      <div><span style="color:var(--muted)">σ_ln</span><br><strong>${sigma_ln.toFixed(3)}</strong></div>
      <div><span style="color:var(--muted)">µ_ln</span><br><strong>${mu_ln.toFixed(3)}</strong></div>
      <div style="grid-column:3/-1"><span style="color:var(--muted)">shape</span><br>
        <strong style="color:var(--accent)">${shape}</strong></div>
    </div>
    <div style="font-size:10px;color:var(--text);line-height:1.5;margin-bottom:4px">${fit}</div>
    ${warning ? '<div style="font-size:10px;color:var(--danger)">'+warning+'</div>' : ''}`;
  desc.style.display = '';
}

function toggleLogScale(canvasId) {
  const ch = window[canvasId+'_chart'];
  if (!ch) return;
  const isLog = window[canvasId+'_logScale'] = !window[canvasId+'_logScale'];
  ch.options.scales.x.type = isLog ? 'logarithmic' : 'linear';
  ch.options.scales.x.ticks.callback = isLog
    ? v => Number(v.toPrecision(2))
    : v => v>=1000 ? (v/1000).toFixed(1)+'k' : v.toFixed(1);
  const btn = document.getElementById(canvasId+'_toggle');
  if (btn) btn.textContent = isLog ? 'Linear x' : 'Log x';
  ch.update();
}

// ── Ca mode switching ────────────────────────────────────────────────
function setCaMode(mode) {
  document.getElementById('ca-manual-block').style.display = mode==='manual' ? '' : 'none';
  document.getElementById('ca-csv-block').style.display    = mode==='csv'    ? '' : 'none';
  document.getElementById('ca-mode-manual').classList.toggle('active', mode==='manual');
  document.getElementById('ca-mode-csv').classList.toggle('active', mode==='csv');
  if (mode === 'manual') fitCaPriorFromManual();
}

function fitCaPriorFromManual() {
  const mean_raw = parseFloat(document.getElementById('ca_conc')?.value);
  const sd_raw   = parseFloat(document.getElementById('ca_conc_sd')?.value);
  const unit     = document.getElementById('ca_unit')?.value || 'ppm';
  const to_ppb   = UNIT_FACTOR[unit] || 1;
  const el = document.getElementById('ca-prior-summary');
  if (!isFinite(mean_raw) || mean_raw <= 0 || !isFinite(sd_raw) || sd_raw <= 0) {
    concentrationPriors['ca'] = null;
    if (el) el.style.display = 'none'; return;
  }
  concentrationPriors['ca'] = fitLognormalFromMeanSd(mean_raw*to_ppb, sd_raw*to_ppb);
  if (el) { el.innerHTML = priorSummaryHTML(concentrationPriors['ca'], 'Ca_aq'); el.style.display = ''; }
}

function fitCaPriorFromCsv() {
  const col  = document.getElementById('ca_csv_col')?.value;
  const unit = document.getElementById('ca_csv_unit')?.value || 'ppm';
  if (!col || !window.caAqCsvData) return;
  const to_ppb = UNIT_FACTOR[unit] || 1;
  const vals = window.caAqCsvData[col].map(v=>parseFloat(v)*to_ppb).filter(v=>isFinite(v)&&v>0);
  concentrationPriors['ca'] = fitLognormal(vals);
  if (!concentrationPriors['ca']) return;
  const mean_ppb = Math.exp(concentrationPriors['ca'].mu_ln + 0.5*concentrationPriors['ca'].sigma_ln**2);
  const caUnit = document.getElementById('ca_unit')?.value || 'ppm';
  document.getElementById('ca_conc').value = (mean_ppb / (UNIT_FACTOR[caUnit]||1)).toFixed(2);
  const el = document.getElementById('ca-prior-summary');
  if (el) { el.innerHTML = priorSummaryHTML(concentrationPriors['ca'],
    `Ca_aq (${vals.length} obs.)`); el.style.display = ''; }
  document.getElementById('ca-dist-chart-wrap').style.display = '';
  renderConcentrationPlot(vals, concentrationPriors['ca'], 'ca-dist-chart', 'ca-dist-desc', 'Ca_aq');
  updateCaHint(); updateAllParamHints();
}

// ── TE_aq mode switching ─────────────────────────────────────────────
function setAqMode(p, mode) {
  const mb = document.getElementById(`${p}-aq-manual-block`);
  const cb = document.getElementById(`${p}-aq-csv-block`);
  if (mb) mb.style.display = mode==='manual' ? '' : 'none';
  if (cb) cb.style.display = mode==='csv'    ? '' : 'none';
  document.getElementById(`${p}-aq-mode-manual`)?.classList.toggle('active', mode==='manual');
  document.getElementById(`${p}-aq-mode-csv`)?.classList.toggle('active', mode==='csv');
  if (mode === 'manual') fitAqPriorFromManual(p);
}

function fitAqPriorFromManual(p) {
  const mean_raw = parseFloat(document.getElementById(`${p}_aq_conc`)?.value);
  const sd_raw   = parseFloat(document.getElementById(`${p}_aq_conc_sd`)?.value);
  const unit     = document.getElementById(`${p}_aq_unit`)?.value || 'ppb';
  const to_ppb   = UNIT_FACTOR[unit] || 1;
  const el = document.getElementById(`${p}-aq-prior-summary`);
  if (!isFinite(mean_raw) || mean_raw <= 0 || !isFinite(sd_raw) || sd_raw <= 0) {
    concentrationPriors[p] = null;
    if (el) el.style.display='none'; return;
  }
  concentrationPriors[p] = fitLognormalFromMeanSd(mean_raw*to_ppb, sd_raw*to_ppb);
  if (el) { el.innerHTML = priorSummaryHTML(concentrationPriors[p],
    `[${p.toUpperCase()}_aq]`); el.style.display = ''; }
}

function fitAqPriorFromCsv(p) {
  const col  = document.getElementById(`${p}_aq_csv_col`)?.value;
  const unit = document.getElementById(`${p}_aq_csv_unit`)?.value || 'ppb';
  if (!col || !window[`${p}AqCsvData`]) return;
  const to_ppb = UNIT_FACTOR[unit] || 1;
  const vals = window[`${p}AqCsvData`][col].map(v=>parseFloat(v)*to_ppb).filter(v=>isFinite(v)&&v>0);
  concentrationPriors[p] = fitLognormal(vals);
  if (!concentrationPriors[p]) return;
  const el = document.getElementById(`${p}-aq-prior-summary`);
  if (el) { el.innerHTML = priorSummaryHTML(concentrationPriors[p],
    `[${p.toUpperCase()}_aq] (${vals.length} obs.)`); el.style.display = ''; }
  const wrap = document.getElementById(`${p}-dist-chart-wrap`);
  if (wrap) wrap.style.display = '';
  const elemLabel = document.getElementById(`${p}_elem`)?.value || p.toUpperCase();
  renderConcentrationPlot(vals, concentrationPriors[p],
    `${p}-dist-chart`, `${p}-dist-desc`, `[${elemLabel}_aq]`);
}

// ── Kd mode toggle: drip rate vs manual ──────────────────────────────────
function setKdMode(prefix, mode) {
  const input  = document.getElementById(prefix + '_Kd_mn');
  const drWrap = document.getElementById(prefix + '_dr_wrap');
  const btnDr  = document.getElementById(prefix + '_kd_mode_dr');
  const btnLit = document.getElementById(prefix + '_kd_mode_lit');
  if (!input) return;

  // Store mode on the input element for other functions to read
  input.dataset.kdMode = mode;

  if (mode === 'driprate') {
    input.readOnly = true;
    input.style.opacity = '0.6';
    input.style.cursor = 'not-allowed';
    if (drWrap) drWrap.style.display = '';
    if (btnDr)  { btnDr.className  = 'btn btn-primary'; btnDr.style.flex  = '1'; }
    if (btnLit) { btnLit.className = 'btn btn-ghost';   btnLit.style.flex = '1'; }
    autoCalcKd(prefix);
  } else {
    input.readOnly = false;
    input.style.opacity = '1';
    input.style.cursor = '';
    if (drWrap) drWrap.style.display = 'none';
    if (btnDr)  { btnDr.className  = 'btn btn-ghost';   btnDr.style.flex  = '1'; }
    if (btnLit) { btnLit.className = 'btn btn-primary';  btnLit.style.flex = '1'; }
  }
  updateParamHints(prefix);
}

function getKdMode(prefix) {
  const input = document.getElementById(prefix + '_Kd_mn');
  return input?.dataset?.kdMode || 'driprate';
}

function getVRef() {
  // Global monitored drip rate for all TEs
  const el = document.getElementById('global_drip_rate');
  return el ? (parseFloat(el.value) || 10) : 10;
}

function autoCalcKd(prefix) {
  if (getKdMode(prefix) !== 'driprate') return;
  const hint = document.getElementById(prefix + '_dr_hint');
  const setHint = msg => { if (hint) hint.innerHTML = msg; };

  const idx    = parseInt(prefix.replace(/[^0-9]/g, '')) - 1;
  const calPct = getCalPct(prefix);
  const obs    = getObsMedianNative(idx, calPct);
  if (obs === null || obs <= 0) {
    setHint('<span style="color:var(--muted)">Upload TE data to auto-calculate Kd</span>');
    return;
  }

  const xf     = parseFloat(document.getElementById(prefix + '_F')?.value) || 0;
  const xl     = parseFloat(document.getElementById(prefix + '_labile')?.value) || 0;
  const keVal  = parseFloat(document.getElementById(prefix + '_K_e')?.value) || 1;
  const kpVal  = parseFloat(document.getElementById(prefix + '_Kp')?.value);
  const kp     = (kpVal === -1 || !isFinite(kpVal))
                 ? (THEO_KP[document.getElementById(prefix+'_elem')?.value] || 1)
                 : kpVal;
  const aqVal  = parseFloat(document.getElementById(prefix + '_aq_conc')?.value);
  const aqUnit = document.getElementById(prefix + '_aq_unit')?.value || 'ppb';
  const aqPPB  = isFinite(aqVal) ? aqVal * (UNIT_FACTOR[aqUnit] || 1) : 0;
  const caVal  = parseFloat(document.getElementById('ca_conc')?.value) || 0;
  const caUnit = document.getElementById('ca_unit')?.value || 'ppb';
  const caPPB  = caVal * (UNIT_FACTOR[caUnit] || 1);
  const vRef   = getVRef();
  const tau    = 60.0 / vRef;

  if (aqPPB <= 0 || caPPB <= 0 || kp <= 0) {
    setHint('<span style="color:#ffa032">Enter valid aq. conc., Ca, and Kp</span>');
    return;
  }

  const nS     = 1.0 - xf;
  const K0_ppm = kp * (xf + xl) * (aqPPB / caPPB) * CA_CALCITE_PPM * keVal;
  if (K0_ppm <= 0) {
    setHint('<span style="color:#ffa032">K₀ ≤ 0 — check Kp and fractions</span>');
    return;
  }

  const dataUnit = teRowData[idx]?.unit || 'ppm';
  const obsPpm   = obs * ((UNIT_FACTOR[dataUnit] || 1) / 1000);
  const ratio    = obsPpm / K0_ppm;
  const fastFloor = K0_ppm * xf;

  if (ratio >= 1.0) {
    setHint(`<span style="color:#ffa032">⚠ Obs (${obsPpm.toFixed(3)} ppm) exceeds K₀ equilibrium (${K0_ppm.toFixed(1)} ppm) — reduce aq. conc or increase Ca</span>`);
    return;
  }

  let E1 = (1 - ratio) / nS;
  let boundaryNote = '';
  if (E1 >= 1.0) {
    E1 = 0.999;
    boundaryNote = ` ⚠ obs near fast-only limit (${fastFloor.toFixed(3)} ppm)`;
  }
  if (E1 <= 0.001) {
    E1 = 0.001;
    boundaryNote = ' ⚠ near equilibrium limit';
  }

  const expKmu = -Math.log(E1) / tau;
  if (expKmu <= 0 || !isFinite(expKmu)) {
    setHint('<span style="color:#ffa032">Cannot solve for Kd at this drip rate</span>');
    return;
  }

  const kdMn = Math.log(expKmu);
  document.getElementById(prefix + '_Kd_mn').value = kdMn.toFixed(3);
  setHint(`τ=${tau.toFixed(1)}s · K₀=${K0_ppm.toFixed(1)} ppm · obs/K₀=${(ratio*100).toFixed(1)}% · ln(Kd)=${kdMn.toFixed(3)}`
    + `<span style="color:#ffa032">${boundaryNote}</span>`);
}


function updateParamHints(prefix) {
  // Auto-calculate Kd first if enabled (before reading Kd_mn for hints)
  autoCalcKd(prefix);
  const idx    = parseInt(prefix.replace(/[^0-9]/g,'')) - 1;
  const kdMn   = parseFloat(document.getElementById(prefix+'_Kd_mn')?.value);
  const kdSd   = parseFloat(document.getElementById(prefix+'_Kd_sd')?.value);
  const xf     = parseFloat(document.getElementById(prefix+'_F')?.value) || 0;
  const xl     = parseFloat(document.getElementById(prefix+'_labile')?.value);
  const keVal  = parseFloat(document.getElementById(prefix+'_K_e')?.value) || 1;
  const aqVal  = parseFloat(document.getElementById(prefix+'_aq_conc')?.value);
  const aqUnit = document.getElementById(prefix+'_aq_unit')?.value||'ppb';
  const aqPPB  = isFinite(aqVal) ? aqVal*(UNIT_FACTOR[aqUnit]||1) : null;
  const molWt  = parseFloat(document.getElementById(prefix+'_mol_wt')?.value) || 58.693;
  const kpVal  = parseFloat(document.getElementById(prefix+'_Kp')?.value);
  // If Kp=-1 use the known lookup value; otherwise use entered value
  const kp     = (kpVal === -1 || !isFinite(kpVal))
                 ? (typeof THEO_KP !== 'undefined'
                    ? (THEO_KP[document.getElementById(prefix+'_elem')?.value] || 1)
                    : 1)
                 : kpVal;
  const caVal  = parseFloat(document.getElementById('ca_conc')?.value) || 0;
  const caUnit = document.getElementById('ca_unit')?.value||'ppb';
  const caPPB  = caVal * (UNIT_FACTOR[caUnit]||1);  // µg/L
  const kd     = isFinite(kdMn) ? Math.exp(kdMn) : null;
  const kdLo   = (isFinite(kdMn)&&isFinite(kdSd)) ? Math.exp(kdMn-kdSd) : null;
  const kdHi   = (isFinite(kdMn)&&isFinite(kdSd)) ? Math.exp(kdMn+kdSd) : null;
  const calPct = getCalPct(prefix);
  const obs    = getObsMedianNative(idx, calPct);

  // Per-TE drip rate → residence time
  const vRef   = getVRef();
  const kdMode = getKdMode(prefix);
  const tau    = 60.0 / vRef;                        // seconds per drip

  // Calibration window hint
  const calHint = document.getElementById(prefix + '_cal_hint');
  if (calHint && kdMode === 'driprate') {
    const row = teRowData[idx];
    if (row && teRawData[row.col] && teRawDepth.length) {
      const proxyVals = (teRawData[row.col]||[]).map(v=>parseFloat(v));
      const finD = teRawDepth.filter(isFinite);
      const dMin = finD.reduce((a,b) => a < b ? a : b, Infinity);
      const dMax = finD.reduce((a,b) => a > b ? a : b, -Infinity);
      const dThr = dMin + (dMax - dMin) * (calPct / 100);
      let nWin = 0;
      for (let i = 0; i < proxyVals.length; i++) {
        if (isFinite(proxyVals[i]) && isFinite(teRawDepth[i]) && teRawDepth[i] <= dThr) nWin++;
      }
      const obsStr = obs !== null ? obs.toFixed(3) : '—';
      calHint.innerHTML = `Using <strong>${nWin}</strong> points (depth ≤ ${dThr.toFixed(1)} ${teDepthUnit}) · `
        + `median = <strong>${obsStr}</strong> ${row.unit || 'ppm'}`;
    }
  }

  // Model physics (matching model.py dr_pdfseries):
  //   K_0 = Kp × (XF+XL) × (aq/ca) × 400000 × K_e     [equilibrium, ppm]
  //   h(V) = K_0 × (1 − nS × E1)                        [at drip rate V]
  //   E1 ≈ exp(−exp(Kd_mn) × τ)                          [narrow-Gaussian approx]
  //   nS = 1 − XF
  // Kp = partition coefficient; Kd = exp(Kd_mn) = OMC dissociation rate (s⁻¹)
  const nS     = 1.0 - xf;
  const K0_ppm = kp * (xf + xl) * (aqPPB > 0 && caPPB > 0 ? aqPPB / caPPB : 0) * CA_CALCITE_PPM * keVal;
  const E1     = kd ? Math.exp(-kd * tau) : null;     // kd = exp(Kd_mn)
  const predEq = K0_ppm;                               // V → 0
  const predKin = (E1 !== null) ? K0_ppm * (1 - nS * E1) : null; // at V_ref

  const h = id => document.getElementById(id);

  // Update per-TE drip rate hint
  const rdh = h(prefix+'_dr_hint');
  if (rdh && kd) {
    const fDiss = 1 - (E1 || 0);
    rdh.textContent = `τ = ${tau.toFixed(2)}s · E₁ = ${(E1||0).toFixed(4)} · ${(fDiss*100).toFixed(1)}% of XL dissociates`;
  }

  // Kd_mn: show Kd value and range
  const kh = h(prefix+'_Kd_mn_hint');
  if (kh) {
    let parts = [];
    if (kd) parts.push(`Kd = ${kd.toExponential(2)} s⁻¹`);
    if (kdLo&&kdHi) parts.push(`1σ: ${kdLo.toExponential(2)}–${kdHi.toExponential(2)}`);
    // In manual mode, show implied Kd_mn from data via kinetic inversion
    if (kdMode === 'manual' && obs && K0_ppm > 0 && nS > 0) {
      const obsPpm = obs * ((UNIT_FACTOR[teRowData[idx]?.unit||'ppm']||1) / 1000);
      const _ratio = obsPpm / K0_ppm;
      const _E1 = (1 - _ratio) / nS;
      if (_E1 > 0 && _E1 < 1) {
        const _impKdMn = Math.log(-Math.log(_E1) / tau);
        parts.push(`implied ln(Kd) @ ${vRef} min⁻¹: ${_impKdMn.toFixed(3)}`);
      }
    }
    kh.textContent = parts.join(' · ');
  }

  // Kd_sd: spread note
  const sh = h(prefix+'_Kd_sd_hint');
  if (sh&&isFinite(kdSd))
    sh.textContent = `±1σ spans ${Math.exp(kdSd).toFixed(2)}× in Kd. Lindeman 2022 default: 1.385.`;

  // aq_conc: show ppb equivalent
  const ah = h(prefix+'_aq_hint');
  if (ah) {
    let parts = [];
    if (aqPPB>0) parts.push(`= ${fmtPPB(aqPPB)} (in ppb)`);
    ah.textContent = parts.join(' · ');
  }

  // Predicted vs observed — using correct model physics
  //   Equilibrium:  K_0 = Kp × (XF+XL) × (aq/ca) × 400000 × K_e
  //   Kinetic:      K_0 × (1 − nS × exp(−Kd×τ))
  const po = h(prefix+'_pred_obs');
  if (po&&analysisMode==='full') {
    const canPredict = K0_ppm > 0 && predKin !== null;
    if (canPredict) {
      const dataUnitNative = teRowData[idx]?.unit || 'ppm';
      const obsPpm = (obs !== null)
        ? obs * ((UNIT_FACTOR[dataUnitNative]||1) / 1000)
        : null;

      let txt = `<div style="display:grid;grid-template-columns:auto 1fr;gap:2px 10px;align-items:baseline">`;
      txt += `<span style="color:var(--muted)">K₀ equilibrium (V→0):</span>`;
      txt += `<span>${predEq.toFixed(3)} ppm</span>`;
      txt += `<span style="color:var(--accent)">Kinetic @ ${vRef} min⁻¹ (τ=${tau.toFixed(1)}s):</span>`;
      txt += `<strong>${predKin.toFixed(3)} ppm</strong>`;
      txt += `<span style="color:var(--muted)">Fast-only limit (V→∞):</span>`;
      txt += `<span>${(K0_ppm * xf).toFixed(4)} ppm</span>`;

      if (obsPpm !== null) {
        const ratioKin = obsPpm / predKin;
        const col  = (ratioKin>5||ratioKin<0.2) ? '#ffa032' : 'var(--accent)';
        txt += `<span style="color:var(--muted)">Observed median:</span>`;
        txt += `<strong>${(obs||0).toFixed(3)} ${dataUnitNative}</strong>`;
        txt += `<span style="color:var(--muted)">Ratio (kinetic):</span>`;
        txt += `<strong style="color:${col}">${ratioKin.toFixed(2)}×</strong>`;
        if (obsPpm > predEq)
          txt += `<span style="grid-column:1/-1;color:#ffa032;margin-top:4px">⚠ Observed exceeds K₀ equilibrium — check Kp, aq_conc, or data units.</span>`;
        else if (obsPpm < K0_ppm * xf)
          txt += `<span style="grid-column:1/-1;color:#ffa032;margin-top:4px">⚠ Observed below fast-only limit — check XF, aq_conc, or data units.</span>`;
        if (ratioKin>5||ratioKin<0.2)
          txt += `<span style="grid-column:1/-1;color:#ffa032">⚠ check units / concentrations</span>`;
      } else {
        txt += `<span style="grid-column:1/-1;color:var(--muted)">Upload TE data to compare</span>`;
      }
      txt += `</div>`;
      po.innerHTML = txt; po.style.display = '';
    } else { po.style.display = 'none'; }
  } else if (po) { po.style.display = 'none'; }

  // Update h(V) diagnostic chart
  renderHvChart(prefix);
}

// ── h(V) diagnostic chart — forward model curve vs data distribution ─────
const _hvCharts = {};  // cache per prefix

// ── Auto-set VMAX from h(V) vs data intersection ────────────────────────
function autoSetVmax() {
  // For the first configured TE, find where h(V) crosses the data p05
  // and set VMAX to 3× that crossing point (or 2× current VMAX if no crossing)
  const hint = document.getElementById('v_max_hint');
  const idx = 0;  // first TE
  const prefix = 'te1';
  const row = teRowData[idx];
  if (!row || !teRawData[row.col]) {
    if (hint) hint.innerHTML = '<span style="color:#ffa032">Upload TE data first</span>';
    return;
  }

  const kdMn  = parseFloat(document.getElementById(prefix+'_Kd_mn')?.value);
  const xf    = parseFloat(document.getElementById(prefix+'_F')?.value) || 0;
  const xl    = parseFloat(document.getElementById(prefix+'_labile')?.value) || 0;
  const keVal = parseFloat(document.getElementById(prefix+'_K_e')?.value) || 1;
  const kpVal = parseFloat(document.getElementById(prefix+'_Kp')?.value);
  const kp    = (kpVal === -1 || !isFinite(kpVal))
                ? (THEO_KP[document.getElementById(prefix+'_elem')?.value] || 1) : kpVal;
  const aqVal = parseFloat(document.getElementById(prefix+'_aq_conc')?.value);
  const aqUnit = document.getElementById(prefix+'_aq_unit')?.value || 'ppb';
  const aqPPB = isFinite(aqVal) ? aqVal * (UNIT_FACTOR[aqUnit] || 1) : 0;
  const caVal = parseFloat(document.getElementById('ca_conc')?.value) || 0;
  const caUnit = document.getElementById('ca_unit')?.value || 'ppb';
  const caPPB = caVal * (UNIT_FACTOR[caUnit] || 1);

  if (!isFinite(kdMn) || aqPPB <= 0 || caPPB <= 0 || kp <= 0) {
    if (hint) hint.innerHTML = '<span style="color:#ffa032">Set model parameters first</span>';
    return;
  }

  const nS = 1.0 - xf;
  const K0 = kp * (xf + xl) * (aqPPB / caPPB) * CA_CALCITE_PPM * keVal;
  const kd = Math.exp(kdMn);

  // Get data p05 (lowest 5th percentile — the smallest concentration the model needs to reach)
  const dataUnit = row.unit || 'ppm';
  const unitScale = (UNIT_FACTOR[dataUnit] || 1) / 1000;
  const vals = (teRawData[row.col] || []).map(v => parseFloat(v) * unitScale).filter(isFinite).sort((a,b) => a-b);
  if (!vals.length) {
    if (hint) hint.innerHTML = '<span style="color:#ffa032">No valid data</span>';
    return;
  }
  const p05 = vals[Math.floor(vals.length * 0.05)];
  const p95 = vals[Math.floor(vals.length * 0.95)];

  // Find V where h(V) = p05 (this is the highest V the model needs)
  // h(V) = K0*(1 - nS*exp(-kd*60/V))
  // Binary search for V where h(V) ≈ p05
  let lo = 0.01, hi = 100000;
  for (let iter = 0; iter < 50; iter++) {
    const mid = (lo + hi) / 2;
    const tau = 60.0 / mid;
    const h = K0 * (1 - nS * Math.exp(-kd * tau));
    if (h > p05) lo = mid; else hi = mid;
  }
  const crossingV = (lo + hi) / 2;
  // Also find V where h(V) = p95 (the lowest V the model needs)
  lo = 0.001; hi = 100000;
  for (let iter = 0; iter < 50; iter++) {
    const mid = (lo + hi) / 2;
    const tau = 60.0 / mid;
    const h = K0 * (1 - nS * Math.exp(-kd * tau));
    if (h > p95) lo = mid; else hi = mid;
  }
  const loV = (lo + hi) / 2;

  // Set VMAX to 3× the crossing point, minimum 2× the data-relevant range
  const suggestedVmax = Math.max(
    Math.ceil(crossingV * 3),
    Math.ceil(loV * 10),
    5  // absolute minimum
  );
  // Round to nice number
  const niceVmax = suggestedVmax <= 10 ? Math.ceil(suggestedVmax)
    : suggestedVmax <= 100 ? Math.ceil(suggestedVmax / 5) * 5
    : Math.ceil(suggestedVmax / 50) * 50;

  document.getElementById('v_max').value = niceVmax;
  if (hint) hint.innerHTML = `h(V) crosses data at V≈${crossingV.toFixed(1)} min⁻¹ → V<sub>max</sub>=${niceVmax}`;
  // Refresh h(V) charts for all TEs
  updateAllParamHints();
}

function renderHvChart(prefix) {
  const canvas = document.getElementById(prefix + '_hv_chart');
  if (!canvas) return;
  const idx = parseInt(prefix.replace(/[^0-9]/g, '')) - 1;

  // Read current parameters
  const kdMn  = parseFloat(document.getElementById(prefix+'_Kd_mn')?.value);
  const xf    = parseFloat(document.getElementById(prefix+'_F')?.value) || 0;
  const xl    = parseFloat(document.getElementById(prefix+'_labile')?.value) || 0;
  const keVal = parseFloat(document.getElementById(prefix+'_K_e')?.value) || 1;
  const kpVal = parseFloat(document.getElementById(prefix+'_Kp')?.value);
  const kp    = (kpVal === -1 || !isFinite(kpVal))
                ? (THEO_KP[document.getElementById(prefix+'_elem')?.value] || 1) : kpVal;
  const aqVal = parseFloat(document.getElementById(prefix+'_aq_conc')?.value);
  const aqUnit = document.getElementById(prefix+'_aq_unit')?.value || 'ppb';
  const aqPPB = isFinite(aqVal) ? aqVal * (UNIT_FACTOR[aqUnit] || 1) : 0;
  const caVal = parseFloat(document.getElementById('ca_conc')?.value) || 0;
  const caUnit = document.getElementById('ca_unit')?.value || 'ppb';
  const caPPB = caVal * (UNIT_FACTOR[caUnit] || 1);

  if (!isFinite(kdMn) || aqPPB <= 0 || caPPB <= 0 || kp <= 0) {
    if (_hvCharts[prefix]) { _hvCharts[prefix].destroy(); delete _hvCharts[prefix]; }
    return;
  }

  const nS = 1.0 - xf;
  const K0 = kp * (xf + xl) * (aqPPB / caPPB) * CA_CALCITE_PPM * keVal;
  const kd = Math.exp(kdMn);
  const vMax = parseFloat(document.getElementById('v_max')?.value) || 100;

  // Compute h(V) for 100 V values
  const N = 100;
  const hvData = [];
  for (let i = 0; i < N; i++) {
    const V = 0.5 + i * (vMax - 0.5) / (N - 1);  // drips/min
    const tau = 60.0 / V;
    const E1 = Math.exp(-kd * tau);
    const h = K0 * (1 - nS * E1);  // ppm
    hvData.push({x: V, y: h});
  }

  // Data distribution from uploaded TE data
  const row = teRowData[idx];
  const dataUnit = row?.unit || 'ppm';
  const unitScale = (UNIT_FACTOR[dataUnit] || 1) / 1000;  // to ppm
  let p25 = null, p50 = null, p75 = null, dMin = null, dMax = null;
  if (row && teRawData[row.col]) {
    const vals = (teRawData[row.col] || []).map(v => parseFloat(v) * unitScale)
                  .filter(isFinite).sort((a,b) => a - b);
    if (vals.length >= 3) {
      p25  = vals[Math.floor(vals.length * 0.25)];
      p50  = vals[Math.floor(vals.length * 0.50)];
      p75  = vals[Math.floor(vals.length * 0.75)];
      dMin = vals[Math.floor(vals.length * 0.05)];
      dMax = vals[Math.floor(vals.length * 0.95)];
    }
  }

  // Build datasets
  const datasets = [
    { label: 'h(V) — forward model',
      data: hvData, type: 'line',
      borderColor: 'rgba(76,201,160,0.9)', borderWidth: 2,
      pointRadius: 0, tension: 0.3, fill: false, yAxisID: 'y' },
  ];

  // Data distribution bands (horizontal lines across full V range)
  if (p50 !== null) {
    // IQR band
    const bandPts = [{x: 0.5, y: p25}, {x: vMax, y: p25}, {x: vMax, y: p75}, {x: 0.5, y: p75}];
    datasets.push({
      label: 'Data IQR (p25–p75)',
      data: [{x: 0.5, y: p25}, {x: vMax, y: p25}],
      type: 'line', borderColor: 'rgba(247,164,64,0.3)', borderWidth: 0,
      pointRadius: 0, fill: { target: {value: p75}, above: 'rgba(247,164,64,0.12)', below: 'rgba(247,164,64,0.12)' },
      yAxisID: 'y',
    });
    // p5-p95 band
    datasets.push({
      label: 'Data range (p5–p95)',
      data: [{x: 0.5, y: dMin}, {x: vMax, y: dMin}],
      type: 'line', borderColor: 'rgba(247,164,64,0.15)', borderWidth: 0,
      pointRadius: 0, fill: { target: {value: dMax}, above: 'rgba(247,164,64,0.06)', below: 'rgba(247,164,64,0.06)' },
      yAxisID: 'y',
    });
    // Median line
    datasets.push({
      label: `Observed median (${p50.toFixed(3)} ppm)`,
      data: [{x: 0.5, y: p50}, {x: vMax, y: p50}],
      type: 'line', borderColor: '#f7a440', borderWidth: 1.5,
      borderDash: [4,3], pointRadius: 0, fill: false, yAxisID: 'y',
    });
  }

  // Equilibrium ceiling line
  datasets.push({
    label: `K\u2080 equilibrium (${K0.toFixed(1)} ppm)`,
    data: [{x: 0.5, y: K0}, {x: vMax, y: K0}],
    type: 'line', borderColor: 'rgba(76,201,160,0.3)', borderWidth: 1,
    borderDash: [6,3], pointRadius: 0, fill: false, yAxisID: 'y',
  });

  // Autoscale y-axis to focus on the data range + relevant h(V) overlap
  // Priority: show the data distribution clearly; K0 ceiling may be off-chart
  let yMin = 0, yMax;
  if (p50 !== null) {
    // Scale to data range with generous padding
    const dataTop = dMax * 1.5;
    const dataBot = Math.max(0, dMin * 0.5);
    // Include h(V) values that are within 3× the data range
    const relevantH = hvData.map(p => p.y).filter(y => y <= dataTop * 2);
    const hMax = relevantH.length ? Math.max(...relevantH) : dataTop;
    yMax = Math.max(dataTop, Math.min(hMax * 1.1, K0 * 1.1));
    yMin = dataBot;
  } else {
    // No data uploaded — show full model range
    yMax = K0 * 1.1;
  }

  if (_hvCharts[prefix]) _hvCharts[prefix].destroy();
  _hvCharts[prefix] = new Chart(canvas, {
    data: { datasets },
    options: {
      animation: false, responsive: true, maintainAspectRatio: false,
      plugins: {
        legend: { labels: { color: '#cdd9e5', font: {size: 9}, boxWidth: 10,
          filter: item => !item.text.startsWith('_') } },
      },
      scales: {
        x: { type: 'linear', min: 0, max: vMax,
             title: { display: true, text: 'Drip rate (min\u207b\u00b9)', color: '#8fa4b5', font: {size: 10} },
             ticks: { color: '#8fa4b5', font: {size: 9} },
             grid: { color: 'rgba(255,255,255,0.04)' } },
        y: { type: 'linear', min: yMin, max: yMax,
             title: { display: true, text: '[TE] calcite (ppm)', color: '#8fa4b5', font: {size: 10} },
             ticks: { color: '#8fa4b5', font: {size: 9} },
             grid: { color: 'rgba(255,255,255,0.04)' } },
      },
    },
  });
}

function updateCaHint() {
  const val  = parseFloat(document.getElementById('ca_conc')?.value);
  const unit = document.getElementById('ca_unit')?.value||'ppb';
  const hint = document.getElementById('ca_hint');
  if (!hint||!isFinite(val)) return;
  const mgL  = val*(UNIT_FACTOR[unit]||1)/1000;
  const flag = mgL>200 ? ' ⚠ unusually high' : mgL<5 ? ' ⚠ unusually low' : '';
  const col  = flag ? '#ffa032' : 'var(--muted)';
  hint.innerHTML = `<span style="color:${col}">${mgL.toFixed(1)} mg/L${flag}</span>`
    +' &nbsp;·&nbsp; typical limestone cave: 20–100 mg/L';
}

function updateAllParamHints() {
  for (let i=1; i<=teRowData.length; i++) updateParamHints(`te${i}`);
  updateCaHint();
}

// Elements that support theoretical Kp (have lattice strain params in model.py kp_theory)
const KP_THEORY_SUPPORTED = ['Cu', 'Ni', 'Co'];

// Literature Kp guidance for elements without Wang & Xu support
// Sources: Huang & Fairchild 2001, Day & Henderson 2013, Tremaine & Froelich 2013
const KP_LITERATURE = {
  'Zn': { range: '4–20',   suggest: 8,   ref: 'Marchitto 2006; Sinclair 2006' },
  'Cd': { range: '1–50',   suggest: 10,  ref: 'Lorens 1981; Rimstidt+ 1998' },
  'Pb': { range: '1–5',    suggest: 2.5, ref: 'Rimstidt+ 1998; Gascoyne 1983' },
  'V':  { range: '0.1–2',  suggest: 0.5, ref: 'Huang+ 2001 (est.)' },
  'Mn': { range: '5–30',   suggest: 10,  ref: 'Lorens 1981; Dromgoole & Walter 1990', caveat: '\u26a0 Mn is strongly redox-controlled — use with caution as drip rate proxy' },
  'Fe': { range: '1–20',   suggest: 4,   ref: 'Dromgoole & Walter 1990 (est.)', caveat: '\u26a0 Fe is strongly redox-controlled and often colloidal — use with caution' },
  'Al': { range: '0.01–1', suggest: 0.1, ref: 'colloid/detrital dominated — Kp poorly defined' },
};

function updateTheoKp(prefix) {
  const kpInput = document.getElementById(prefix + '_Kp');
  const theoField = document.getElementById(prefix + '_Kp_theo');
  const val = parseFloat(kpInput.value);
  const elem = document.getElementById(prefix + '_elem').value;

  if (val === -1 || kpInput.value.trim() === '-1') {
    if (THEO_KP[elem] !== undefined) {
      // Empirical from Lindeman — best available
      theoField.value = 'Kp = ' + THEO_KP[elem] + ' (Lindeman 2022)';
      theoField.style.opacity = '0.7';
      theoField.style.color = '';
    } else if (KP_THEORY_SUPPORTED.includes(elem)) {
      // Wang & Xu lattice strain — will auto-calculate
      theoField.value = 'auto-calc via Wang & Xu (2001)';
      theoField.style.opacity = '0.7';
      theoField.style.color = '';
    } else if (KP_LITERATURE[elem]) {
      // Has literature guidance — show it and auto-fill
      const lit = KP_LITERATURE[elem];
      kpInput.value = lit.suggest;
      theoField.value = 'Range ' + lit.range + ' (' + lit.ref + ')';
      if (lit.caveat) theoField.value = lit.caveat + ' | ' + theoField.value;
      theoField.style.opacity = '0.9';
      theoField.style.color = '#f7a440';
    } else {
      theoField.value = '\u26a0 No Kp data for ' + elem + ' \u2014 enter manually';
      theoField.style.opacity = '1';
      theoField.style.color = '#ffa032';
    }
  } else {
    // User has set a manual value — show literature range if available
    if (KP_LITERATURE[elem]) {
      const lit = KP_LITERATURE[elem];
      theoField.value = (lit.caveat ? lit.caveat + ' | ' : '') + 'lit. range: ' + lit.range + ' (' + lit.ref + ')';
      theoField.style.opacity = '0.6';
      theoField.style.color = '';
    } else if (THEO_KP[elem] !== undefined) {
      theoField.value = 'empirical: ' + THEO_KP[elem] + ' (Lindeman 2022)';
      theoField.style.opacity = '0.5';
      theoField.style.color = '';
    } else {
      theoField.value = '';
      theoField.style.opacity = '0.4';
      theoField.style.color = '';
    }
  }
}

// ── Labile fraction auto-calculation ────────────────────────────────────────────
function updateLabile(prefix) {
  const xF = parseFloat(document.getElementById(prefix + '_F').value) || 0;
  const xI = parseFloat(document.getElementById(prefix + '_InertF').value) || 0;
  const xL = Math.max(0, 1.0 - xI - xF);
  document.getElementById(prefix + '_labile').value = xL.toFixed(4);
  updateParamHints(prefix);
  // Warn if fractions exceed 1
  const warn = (xI + xF) > 1.0;
  document.getElementById(prefix + '_F').style.borderColor    = warn ? 'var(--danger)' : '';
  document.getElementById(prefix + '_InertF').style.borderColor = warn ? 'var(--danger)' : '';
}

// ── Molecular weight lookup on element change ────────────────────────────────
const MOL_WT = {
  'Ni': 58.693, 'Co': 58.933, 'Cu': 63.546, 'V': 50.942,
  'Zn': 65.38,  'Cd': 112.411,'Al': 26.982, 'Pb': 207.2,
  'Mn': 54.938, 'Fe': 55.845,
  'La': 138.905,'Ce': 140.116,'other': null
};
// Speleothem-specific inorganic Kd values from Lindeman et al. GCA 2022
const KP_DEFAULT = {
  'Ni': 1.1, 'Co': 4.4, 'Cu': 44,
  'Zn': 8, 'Cd': 10, 'Pb': 2.5, 'V': 0.5,
  'Mn': 10, 'Fe': 4, 'Al': 0.1, 'other': -1
};
function updateMolWt(prefix, elem) {
  const wt = MOL_WT[elem];
  if (wt !== null && wt !== undefined) {
    document.getElementById(prefix + '_mol_wt').value = wt;
  }
  // Update aqueous concentration label to show current element
  const aqLabel = document.getElementById(prefix + '_aq_elem_label');
  if (aqLabel) aqLabel.textContent = elem;
  // Also update dynamic spans in dropzone labels
  document.querySelectorAll('.' + prefix + '-aq-elem-dyn').forEach(el => el.textContent = elem);
  // Set Kp default for this element
  const kp = KP_DEFAULT[elem] !== undefined ? KP_DEFAULT[elem] : -1;
  document.getElementById(prefix + '_Kp').value = kp;
  updateTheoKp(prefix);
  // Propagate element-specific kinetic defaults (Kd, fractions, aq_conc)
  const d = ELEM_DEFAULTS[elem] || ELEM_DEFAULT_FALLBACK;
  const fields = {Kd_mn: d.Kd_mn, Kd_sd: d.Kd_sd, K_e: d.K_e, F: d.F, InertF: d.InertF, aq_conc: d.aq_conc};
  for (const [f, v] of Object.entries(fields)) {
    const el = document.getElementById(prefix + '_' + f);
    if (el) el.value = v;
  }
  // Set aq_conc_sd if available, otherwise clear it
  const sdEl = document.getElementById(prefix + '_aq_conc_sd');
  if (sdEl) sdEl.value = d.aq_conc_sd || '';
  updateLabile(prefix);
  updateParamHints(prefix);
  // Refit stochastic prior from new defaults
  if (typeof fitAqPriorFromManual === 'function') fitAqPriorFromManual(prefix);
}

// ── Utilities ────────────────────────────────────────────────────────────────
function fmtTime(s) {
  if (s < 60) return Math.round(s) + 's';
  return Math.round(s/60) + 'min';
}
function escHtml(s) {
  return s.replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
}

// ── Initialise on page load ───────────────────────────────────────────────
document.addEventListener('DOMContentLoaded', () => {
  console.log('%c PaleoDripRates v44 loaded ', 'background:#4cc9a0;color:#0d1117;font-weight:bold;padding:2px 8px;border-radius:4px');
  renderTEParamCards();
  setAnalysisMode(analysisMode);
  // Fit stochastic priors from default concentrations
  fitCaPriorFromManual();
  for (let i = 1; i <= teRowData.length; i++) fitAqPriorFromManual('te' + i);
});
</script>
</body>
</html>
'''

@app.route('/')
def index():
    resp = Response(HTML, mimetype='text/html')
    resp.headers['Cache-Control'] = 'no-cache, no-store, must-revalidate'
    resp.headers['Pragma'] = 'no-cache'
    resp.headers['Expires'] = '0'
    return resp


@app.route('/upload', methods=['POST'])
def upload():
    """Accept CSV uploads for depth/age, trace elements, and isotopes."""
    results = {}
    for key in ['depth_age', 'trace_elem1', 'trace_elem2', 'isotope1', 'isotope2',
                'ca_aq', 'te1_aq', 'te2_aq', 'te3_aq', 'te4_aq']:
        f = request.files.get(key)
        if f and f.filename:
            path = os.path.join(UPLOAD_FOLDER, key + '.csv')
            f.save(path)
            try:
                df = pd.read_csv(path)
                result_entry = {
                    'filename': f.filename,
                    'columns': list(df.columns),
                    'rows': len(df),
                    'preview': df.head(3).to_dict(orient='records'),
                }
                # For depth_age, trace_elem1, and monitoring CSVs, return full data
                if key in ('depth_age', 'trace_elem1') or key.endswith('_aq'):
                    result_entry['data'] = df.to_dict(orient='list')
                results[key] = result_entry
            except Exception as e:
                results[key] = {'error': str(e)}
    return jsonify(results)


@app.route('/preprocess', methods=['POST'])
def preprocess():
    """Block-average + sigma-clip the uploaded trace_elem1.csv and overwrite it."""
    try:
        p = request.get_json()
        dep_col  = p.get('dep_col', '')
        target_n = int(p.get('target_n', 500))
        sigma    = float(p.get('sigma', 3.0))
        win_size = int(p.get('win_size', 50))

        path = os.path.join(UPLOAD_FOLDER, 'trace_elem1.csv')
        df   = pd.read_csv(path)
        orig = len(df)

        # 1. Windowed (or global) sigma-clip per proxy column
        num_cols = df.select_dtypes(include='number').columns.tolist()
        n_removed = 0
        for col in num_cols:
            if col == dep_col:
                continue
            s = df[col].astype(float)
            if win_size > 0:
                # Centred rolling window — min_periods=3 avoids NaN-filled edges
                roll = s.rolling(window=win_size, center=True, min_periods=3)
                mu_w = roll.mean()
                sd_w = roll.std()
            else:
                mu_w = pd.Series(np.nanmean(s), index=s.index)
                sd_w = pd.Series(np.nanstd(s),  index=s.index)
            mask = (s - mu_w).abs() > sigma * sd_w
            mask = mask.fillna(False)
            n_removed += int(mask.sum())
            df.loc[mask, col] = np.nan

        # 2. Drop rows where all proxy cols are NaN
        proxy_cols = [c for c in num_cols if c != dep_col]
        df = df.dropna(subset=proxy_cols, how='all').reset_index(drop=True)

        # 3. Block-average to target_n rows
        if target_n < len(df) and target_n >= 10:
            block = max(1, len(df) // target_n)
            groups = df.groupby(df.index // block)
            df = groups.mean(numeric_only=True).reset_index(drop=True)

        df.to_csv(path, index=False)
        return jsonify({'rows': len(df), 'original_rows': orig,
                        'removed_outliers': n_removed})
    except Exception as e:
        return jsonify({'error': str(e)}), 400


@app.route('/run', methods=['POST'])
def run():
    """Kick off the model run in a background thread."""
    global run_state
    if run_state['running']:
        return jsonify({'error': 'A run is already in progress.'}), 400

    params = request.get_json()
    run_state = {
        'running': True, 'progress': 0, 'stage': 'Initialising',
        'log': [], 'error': None, 'done': False, 'outputs': [],
    }
    t = threading.Thread(target=_run_model, args=(params,), daemon=True)
    t.start()
    return jsonify({'status': 'started'})


@app.route('/status')
def status():
    return jsonify({
        'running':  run_state['running'],
        'progress': run_state['progress'],
        'stage':    run_state['stage'],
        'log':      run_state['log'][-50:],   # last 50 lines
        'error':    run_state['error'],
        'done':     run_state['done'],
        'outputs':  run_state['outputs'],
    })


@app.route('/download/<filename>')
def download(filename):
    safe = os.path.basename(filename)
    path = os.path.join(OUTPUT_FOLDER, safe)
    if not os.path.exists(path):
        return jsonify({'error': 'File not found'}), 404
    return send_file(path, as_attachment=True)

@app.route('/data/<filename>')
def serve_data(filename):
    """Serve output files inline (for fetch/AJAX) without Content-Disposition: attachment."""
    safe = os.path.basename(filename)
    path = os.path.join(OUTPUT_FOLDER, safe)
    if not os.path.exists(path):
        return jsonify({'error': 'File not found'}), 404
    if safe.endswith('.json'):
        with open(path) as f:
            return Response(f.read(), mimetype='application/json')
    return send_file(path, as_attachment=False)


# ── Model runner (background thread) ────────────────────────────────────────

def _run_model(params):
    """
    Drives the full pipeline using parameters from the UI.
    Mirrors Drip_rate_parallel.py but reads from uploaded CSVs
    rather than an Excel workbook.
    """
    try:
        # Add the parent directory (where bayprox, model.py etc. live) to path
        parent = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
        if parent not in sys.path:
            sys.path.insert(0, parent)

        import bayprox
        dat  = bayprox.data
        ad   = bayprox.agedepth
        prx  = bayprox.proxyrecord

        from scipy.interpolate import interp1d
        import model
        from params import (BAYPROX_SAMPLING_RES as sampling_res,
                            BAYPROX_MAX_DEGREE as max_degree,
                            BAYPROX_SOLVER_METHOD as method)

        def _unit_to_ppb(unit):
            """Return multiplier to convert data in `unit` to ppb.
            Defined here (outside cache/no-cache branches) so it is always
            available to _aq() and make_te() regardless of whether the proxy
            record is loaded from cache or computed fresh.
            """
            return {'ppb': 1.0, 'ppm': 1000.0,
                    'ug/g': 1000.0, 'mg/kg': 1000.0}.get(str(unit).strip(), 1.0)

        # Ca fraction in CaCO₃ by mass: 40.078/100.087 ≈ 40% = 400,000 ppm
        CA_CALCITE_PPM = 400000

        # ── Check whether we can use the cached proxy record ────────────────
        proxy_pkl = os.path.join(OUTPUT_FOLDER, 'ProxyRecord.pkl')
        use_cache = params.get('use_cached_proxy') and os.path.exists(proxy_pkl)

        if use_cache:
            # ── Fast path: load cached PDist objects, skip data/BayProX ────
            _set_stage('Loading cached proxy record', 15)
            log('Found ProxyRecord.pkl — skipping data loading and BayProX.')
            with open(proxy_pkl, 'rb') as f:
                cached = pickle.load(f)
            # Support both the web app format [TE1, TE2, iso]
            # and the original script format (nested dataAll list)
            # Determine how many items are TEs vs isotope
            # Convention: last item is iso if has_iso was True when saved,
            # otherwise all items are TEs. We infer by checking params.
            _te_list_len = len(params.get('te_list') or [1])  # at least 1
            if isinstance(cached, list) and len(cached) >= 1 and not isinstance(cached[0], list):
                # New web-app format
                PDist_TEs = cached[:_te_list_len]
                PDist_iso = cached[_te_list_len] if len(cached) > _te_list_len else None
                PDist_TE1 = PDist_TEs[0]
                PDist_TE2 = PDist_TEs[1] if len(PDist_TEs) > 1 else None
                log(f'Loaded ProxyRecord.pkl ({len(PDist_TEs)} TE(s), iso={PDist_iso is not None})')
            else:
                # Original script format
                PDist_TEs = [cached[1][0][0][3], cached[1][1][0][3]]
                PDist_TE1, PDist_TE2 = PDist_TEs
                PDist_iso = cached[2][0][0][3]
                log('Loaded original script style ProxyRecord.pkl')

        else:
            # ── Full path: load CSVs, run BayProX ───────────────────────────
            from drip_rate_util import detect_outliers, preremvar_optimize_residual
            from params import (OUTLIER_WINDOW_SIZE as WS,
                                OUTLIER_ZERO_TOL as ZERO_TOL,
                                OUTLIER_PDF_TOL  as PDF_TOL)

            # ── 1. Load data ─────────────────────────────────────────────────
            _set_stage('Loading data', 2)


            def load(key, col_depth, col_proxy, unit='ppb'):  # unit kept for back-compat only
                path = os.path.join(UPLOAD_FOLDER, key + '.csv')
                df = pd.read_csv(path)
                x = df[col_depth].to_numpy(dtype=float)
                # Solid proxy data is NOT unit-converted here — model.dr_pdfseries
                # expects data in the same units as when the model was calibrated.
                # The unit selector is used only for hints and sanity checks.
                y = df[col_proxy].to_numpy(dtype=float)
                mask = ~np.isnan(x) & ~np.isnan(y)
                return x[mask], y[mask]

            depth_age_path = os.path.join(UPLOAD_FOLDER, 'depth_age.csv')
            da_df = pd.read_csv(depth_age_path)
            dating_depth     = da_df[params['col_depth']].to_numpy(dtype=float)
            dating_age       = da_df[params['col_age']].to_numpy(dtype=float)
            dating_age_error = da_df[params['col_age_err']].to_numpy(dtype=float) / 2.

            # Build TE list from te_list param (new dynamic system)
            # Falls back to legacy te1_col_*/te2_col_* if te_list absent
            _te_list = params.get('te_list', [])
            if not _te_list:
                _te_list = [{'col_depth': params.get('te1_col_depth',''),
                              'col_proxy': params.get('te1_col_proxy',''),
                              'unit':      params.get('te1_unit','ppb')}]
                if params.get('te2_col_proxy') and \
                   params.get('te2_col_proxy') != params.get('te1_col_proxy'):
                    _te_list.append({'col_depth': params.get('te2_col_depth',''),
                                     'col_proxy': params.get('te2_col_proxy',''),
                                     'unit':      params.get('te2_unit','ppb')})

            TE_data = []  # list of (x, y, row_key)
            for _i, _te in enumerate(_te_list):
                _rk = f'te{_i+1}'
                _x, _y = load('trace_elem1', _te['col_depth'], _te['col_proxy'],
                               _te.get('unit', 'ppb'))
                TE_data.append((_x, _y, _rk))
                log(f'Loaded {_rk}: {len(_x)} points ({_te["col_proxy"]} @ {_te.get("unit","ppb")})' )

            # convenience aliases used by later stages
            x_TE1, y_TE1 = TE_data[0][0], TE_data[0][1]

            _iso_path = os.path.join(UPLOAD_FOLDER, 'isotope1.csv')
            has_iso = (os.path.isfile(_iso_path)
                       and bool(params.get('iso_col_depth'))
                       and bool(params.get('iso_col_proxy')))
            if has_iso:
                x_iso, y_iso = load('isotope1', params['iso_col_depth'], params['iso_col_proxy'])
            else:
                x_iso, y_iso = None, None
                log('Isotope file not found — running without isotope')

            log(f'Loaded depth/age: {len(dating_depth)} points')
            log(f'Loaded {len(TE_data)} trace element(s)')
            if has_iso: log(f'Loaded isotope: {len(x_iso)} points')

            # ── 1b. Depth unit harmonisation ──────────────────────────────────
            # Ensure TE depths are in the same unit as depth_age depths.
            _da_unit = params.get('da_depth_unit', 'cm')
            _te_unit = params.get('te_depth_unit', 'cm')
            _depth_scale = 1.0
            if _da_unit != _te_unit:
                if _te_unit == 'mm' and _da_unit == 'cm':
                    _depth_scale = 0.1   # mm → cm
                elif _te_unit == 'cm' and _da_unit == 'mm':
                    _depth_scale = 10.0  # cm → mm
                log(f'⚠ Depth unit mismatch: depth_age={_da_unit}, TE={_te_unit} → '
                    f'scaling TE depths by {_depth_scale}')
                TE_data = [(_x * _depth_scale, _y, _rk) for _x, _y, _rk in TE_data]
                x_TE1 = TE_data[0][0]  # refresh alias
                if has_iso and x_iso is not None:
                    x_iso = x_iso * _depth_scale
                    log(f'  Also scaled isotope depths by {_depth_scale}')
            else:
                log(f'Depth units match: {_da_unit} (no conversion needed)')
            log(f'Depth_age range: {dating_depth.min():.2f}–{dating_depth.max():.2f} {_da_unit}')
            log(f'TE1 depth range: {x_TE1.min():.2f}–{x_TE1.max():.2f} {_da_unit} (after conversion)')

            # ── 2. Unified depth grid ────────────────────────────────────────
            _set_stage('Building unified depth grid', 5)
            if has_iso:
                unified_depth = np.sort(np.r_[x_TE1, x_iso])
                xmin = max(x_TE1.min(), x_iso.min())
                xmax = min(x_TE1.max(), x_iso.max())
            else:
                unified_depth = np.sort(x_TE1)
                xmin, xmax = x_TE1.min(), x_TE1.max()
            unified_depth = unified_depth[(unified_depth >= xmin) & (unified_depth <= xmax)]
            log(f'Unified depth: {len(unified_depth)} points ({xmin:.2f} – {xmax:.2f} cm)')

            # ── 3. Outlier removal ───────────────────────────────────────────
            _set_stage('Removing outliers', 8)
            _ui_ws = int(params.get('outlier_win_size', WS))
            def clean(x, y, name):
                y_ = y.copy()
                idx = detect_outliers(name, x, y, winsize=_ui_ws,
                                      zero_tol=ZERO_TOL, pdf_tol=PDF_TOL)
                kk = int(_ui_ws / 2)
                for ii in idx:
                    y_[ii] = np.median(y[max(0,ii-kk):ii+kk])
                log(f'{name}: removed {len(idx)} outliers (window={_ui_ws}pts)')
                return y_

            TE_data = [(_x, clean(_x, _y, _rk.upper()), _rk)
                       for _x, _y, _rk in TE_data]
            y_TE1 = TE_data[0][1]  # keep alias fresh
            if has_iso:
                y_iso = clean(x_iso, y_iso, 'ISO')

            # ── 4. Age model ─────────────────────────────────────────────────
            _set_stage('Building age model', 10)
            CALAGE_MIN = int(params['calage_min'])
            CALAGE_MAX = int(params['calage_max'])
            calage = np.arange(CALAGE_MIN, CALAGE_MAX, sampling_res)

            info = dat.SampleInfo(name=params['station_name'],
                                  datatype='none', archive='Speleothem')
            DT1 = dat.DatingTable(depth=dating_depth, age=dating_age,
                                  ageerror=dating_age_error,
                                  datingmethod='U/Th', sampleinfo=info)
            DT1.calibration()

            # ── 5. BayProX proxy record ──────────────────────────────────────
            _set_stage('Computing proxy record', 15)
            log('Starting BayProX proxy record computation …')

            PD_unif = dat.ProxyDepth(depth=unified_depth,
                                     proxy=np.zeros(len(unified_depth)),
                                     proxyerror=0., sampleinfo=info)
            pr_rem_var = np.linspace(0., 0.20, 100)
            pr_rem_res = np.zeros(100)
            for i in range(100):
                pr_rem_res[i] = preremvar_optimize_residual(
                    pr_rem_var[i], DT1, PD_unif, calage, max_degree, method, -1)
            pre_remvar = pr_rem_var[np.argmin(pr_rem_res)] * DT1.ageerror.std()
            log(f'Optimal prior remainder variance: {pre_remvar:.4f}')
            _set_stage('Computing proxy distributions', 25)

            def get_pdist(x, y, dtype, limits_mode='te'):
                f_y = interp1d(x, y, kind='linear',
                               bounds_error=False, fill_value=np.nan)
                y_u = f_y(unified_depth)
                info_i = dat.SampleInfo(name=params['station_name'],
                                        datatype=dtype, archive='Speleothem')
                PD = dat.ProxyDepth(depth=unified_depth, proxy=y_u,
                                    proxyerror=0., sampleinfo=info_i)
                DWF = ad.DWF(DT1)
                DWF(PD.depth, calBP_axis=calage, max_degree=max_degree,
                    pre_remvar=pre_remvar, method=method, verbose=0)
                PDist = prx.ProxyDistributions(DWF)
                if limits_mode == 'te':
                    _y_valid = PD.proxy[np.isfinite(PD.proxy)]
                    p25, p50, p75 = np.percentile(_y_valid, [25, 50, 75])
                    limits = [0., p50 + 10. * (p75 - p25)]
                    log(f'  PDist limits ({dtype}): p25={p25:.4f} p50={p50:.4f} p75={p75:.4f} '
                        f'→ [{limits[0]:.4f}, {limits[1]:.4f}]  '
                        f'data range: [{np.nanmin(_y_valid):.4f}, {np.nanmax(_y_valid):.4f}]')
                else:
                    limits = PDist.get_limits(DWF, PD, 3.)
                    log(f'  PDist limits ({dtype}): auto → [{limits[0]:.4f}, {limits[1]:.4f}]')
                PDist.get_pdf(DWF, PD, res=500, limits=limits)
                return PDist

            PDist_TEs = []
            _n_te = len(TE_data)
            for _i, (_x, _y, _rk) in enumerate(TE_data):
                _pct_lo = 25 + int(20 * _i / _n_te)
                _pct_hi = 25 + int(20 * (_i+1) / _n_te)
                _set_stage(f'Computing proxy distributions ({_rk.upper()})', _pct_lo)
                log(f'Computing {_rk.upper()} proxy distribution …')
                PDist_TEs.append(get_pdist(_x, _y, 'concentration', 'te'))
            PDist_TE1 = PDist_TEs[0]  # keep alias
            PDist_TE2 = PDist_TEs[1] if len(PDist_TEs) > 1 else None

            if has_iso:
                _set_stage('Computing proxy distributions (isotope)', 60)
                log('Computing isotope proxy distribution …')
                PDist_iso = get_pdist(x_iso, y_iso, 'isotope', 'iso')
            else:
                PDist_iso = None

            _pkl_list = PDist_TEs + ([PDist_iso] if PDist_iso is not None else [])
            with open(proxy_pkl, 'wb') as f:
                pickle.dump(_pkl_list, f)
            log('Proxy record saved to ProxyRecord.pkl')

        # ── 6. Drip rate PDFs ───────────────────────────────────────────────
        _set_stage('Computing drip rate PDF (TE1)', 65)
        log('Computing drip rate PDF for TE1 …')

        from driprates_stochastic import driprates

        # ── Apply user VMAX / VRES overrides ─────────────────────────────
        # model.py and model_stochastic.py do `from params import VMAX, VMIN, VRES`
        # at module load — they hold local copies. Monkeypatch both modules directly.
        import params as _params_mod
        import model_stochastic as _ms_mod
        _user_vmax = float(params.get('v_max', 100))
        _user_vres = int(params.get('v_res', 5000))
        _user_vmin = _user_vmax / _user_vres  # match params.py: VMIN = VMAX / VRES
        for _mod in [_params_mod, model, _ms_mod]:
            _mod.VMAX = _user_vmax
            _mod.VMIN = _user_vmin
            _mod.VRES = _user_vres
        log(f'V grid: VMAX={_user_vmax}, VRES={_user_vres}, VMIN={_user_vmin:.6f}')

        # ── Monkeypatch te_pdfseries to prevent proxyspan in-place mutation ──
        # model.te_pdfseries does X_s = PDist.proxyspan; X_s /= (1E3 * mol_wt)
        # which mutates PDist.proxyspan in place. Without deepcopy (which
        # corrupts BayProX objects), subsequent calls divide again → collapse.
        # Fix: wrap with a lightweight PDist copy that isolates proxyspan.
        _orig_te_pdfseries = model.te_pdfseries

        class _PDist_Wrap:
            """Lightweight proxy for ProxyDistributions with isolated proxyspan."""
            pass

        def _safe_te_pdfseries(TE):
            _w = _PDist_Wrap()
            _w.calbp      = TE['PDist'].calbp
            _w.proxyspan  = TE['PDist'].proxyspan.copy()  # isolated copy
            _w.pdfmat     = TE['PDist'].pdfmat
            _TE2 = dict(TE)
            _TE2['PDist'] = _w
            return _orig_te_pdfseries(_TE2)

        model.te_pdfseries = _safe_te_pdfseries
        log('Monkeypatched model.te_pdfseries (proxyspan copy isolation)')

        # ── Inline mini-test: reproduce dr_pdfseries for 1 timestep ──────
        # This bypasses model.py entirely to pinpoint NaN source.
        def _inline_test(pd_obj, te_dict, kd_mn, kd_sd, k_e):
            from params import KRES, KLIM
            VMIN, VMAX, VRES = _user_vmin, _user_vmax, _user_vres
            _ps = pd_obj.proxyspan.copy()
            _pm = pd_obj.pdfmat
            _mw = te_dict['mol_wt']
            log(f'  [TEST] proxyspan raw: [{_ps[0]:.6g}, {_ps[-1]:.6g}], len={len(_ps)}')
            _ps_mol = _ps / (1E3 * _mw)
            log(f'  [TEST] proxyspan mol/kg: [{_ps_mol[0]:.6g}, {_ps_mol[-1]:.6g}]')
            # Build interp1d for timestep 0
            from scipy.interpolate import interp1d as _i1d
            _f0 = _i1d(_ps_mol, _pm[:, 0], bounds_error=False, fill_value=0.)
            # Aqueous
            _inertF = te_dict['InertF']
            _Xa = (1.0 - _inertF) * te_dict['aq_conc'] / (1E6 * _mw)
            _Ya = te_dict['ca_conc'] / (1E6 * 40.078)
            _Ys = 400. / 40.078
            _Kp = te_dict['Kp']
            _K0 = _Kp * _Ys * (_Xa / _Ya) * k_e
            _nF = te_dict['F']
            _nS = 1. - _nF
            _k_mu, _k_sd = kd_mn, kd_sd
            log(f'  [TEST] Xa={_Xa:.6e} Ya={_Ya:.6e} K0={_K0:.6e} nS={_nS}')
            # V span
            _V = np.linspace(VMIN, VMAX, VRES) / 60.  # drips/sec
            # k span
            _k = np.linspace(_k_mu - KLIM*_k_sd, _k_mu + KLIM*_k_sd, KRES)
            _k_rsw = 0.5 * np.r_[_k[1]-_k[0], _k[2:]-_k[:-2], _k[-1]-_k[-2]]
            log(f'  [TEST] V range: [{_V[0]:.6e}, {_V[-1]:.6e}] drips/s  k range: [{_k[0]:.2f}, {_k[-1]:.2f}]')
            # Compute E_1 and h(v) for 5 test V values
            _test_idxs = [0, VRES//4, VRES//2, 3*VRES//4, VRES-1]
            for _vi in _test_idxs:
                _v = _V[_vi]
                # E_1 = integral of gaussian(k)*exp(-exp(k)/v) dk
                _ik = np.exp(- np.exp(_k) / _v)   # exp(-K_d/v)
                _gk = np.exp(-(_k - _k_mu)**2 / (2*_k_sd**2)) / (_k_sd * np.sqrt(2*np.pi))
                _E1 = (_ik * _gk * _k_rsw).sum()
                _hv = _K0 * (1. - _nS * _E1)
                # E_2 for h'(v)
                _ik2 = _ik * np.exp(_k) / _v**2
                _E2 = (_ik2 * _gk * _k_rsw).sum()
                _hp = -_K0 * _nS * _E2
                _pdf_val = _f0(_hv)
                _raw_pdf = - _pdf_val * _hp
                _vpm = _v * 60.
                log(f'  [TEST] V={_vpm:.1f}/min: E1={_E1:.6f} h={_hv:.6e} hp={_hp:.6e} '
                    f'pdf(h)={_pdf_val:.4f} raw_pdf={_raw_pdf:.6e} '
                    f'{"✓ in range" if 0 < _hv < _ps_mol[-1] else "✗ OUT OF RANGE"}')
            # Full V computation
            _nk, _nv = len(_k), len(_V)
            _ka = _k.repeat(_nv).reshape(_nk, _nv).T
            _va = _V.repeat(_nk).reshape(_nv, _nk)
            _t0 = -np.exp(_ka) / _va
            _t1 = np.exp(_t0)
            _t2 = np.exp(-(_ka - _k_mu)**2 / (2*_k_sd**2)) / (_k_sd*np.sqrt(2*np.pi))
            _fk = _t1 * _t2
            _E1_all = (_fk * _k_rsw).sum(axis=1)
            _h_all = _K0 * (1. - _nS * _E1_all)
            _t3 = np.exp(_ka) / _va**2
            _fk2 = _t1 * _t2 * _t3
            _E2_all = (_fk2 * _k_rsw).sum(axis=1)
            _hp_all = -_K0 * _nS * _E2_all
            log(f'  [TEST] h_all: min={np.nanmin(_h_all):.6e} max={np.nanmax(_h_all):.6e} '
                f'NaN={np.isnan(_h_all).sum()} inf={np.isinf(_h_all).sum()}')
            log(f'  [TEST] hp_all: min={np.nanmin(_hp_all):.6e} max={np.nanmax(_hp_all):.6e} '
                f'NaN={np.isnan(_hp_all).sum()} inf={np.isinf(_hp_all).sum()}')
            _pdf0 = _f0(_h_all)
            log(f'  [TEST] pdf0(h): min={np.nanmin(_pdf0):.6e} max={np.nanmax(_pdf0):.6e} '
                f'NaN={np.isnan(_pdf0).sum()} nonzero={np.count_nonzero(_pdf0)}')
            _raw = -_pdf0 * _hp_all
            _raw_sum = np.nansum(_raw * (0.5 * np.r_[_V[1]-_V[0], _V[2:]-_V[:-2], _V[-1]-_V[-2]]))
            log(f'  [TEST] raw_pdf: min={np.nanmin(_raw):.6e} max={np.nanmax(_raw):.6e} '
                f'NaN={np.isnan(_raw).sum()} sum={_raw_sum:.6e}')
            if _raw_sum > 0:
                log(f'  [TEST] ✓ INLINE COMPUTATION PRODUCES VALID PDF — model.py issue suspected')
            else:
                log(f'  [TEST] ✗ Inline also fails — check K0, proxyspan, or unit mismatch')
        _analysis_mode = params.get('analysis_mode', 'full')
        log(f'Analysis mode: {_analysis_mode}')

        # ── Diagnostic: dump raw TE params as received from UI ────────────
        for _di in range(len(params.get('te_list') or [{'col_proxy':'?'}])):
            _dp = f'te{_di+1}'
            _dkeys = ['elem','mol_wt','Kp','Kd_mn','Kd_sd','K_e','cal_pct','F','InertF','labile','aq_conc','aq_unit']
            _dvals = {k: params.get(f'{_dp}_{k}', '?') for k in _dkeys}
            log(f'Raw params {_dp}: {_dvals}')
        log(f'Raw params ca_conc={params.get("ca_conc","?")} ca_unit={params.get("ca_unit","?")}')

        # Build a lookup of observed median concentrations (ppb) per TE row key
        # Used in semi mode to auto-scale aq_conc so the PDF has support in data range.
        # Load proxy medians (native units, unconverted) for semi-quant auto-scaling
        _obs_median_native = {}
        try:
            _te_csv = os.path.join(UPLOAD_FOLDER, 'trace_elem1.csv')
            if os.path.isfile(_te_csv):
                _te_df = pd.read_csv(_te_csv)
                for _i, _te in enumerate(params.get('te_list', []) or
                                         [{'col_proxy': params.get('te1_col_proxy',''),
                                           'unit': params.get('te1_unit','ppm')}]):
                    _col = _te.get('col_proxy','')
                    if _col and _col in _te_df.columns:
                        _vals = pd.to_numeric(_te_df[_col], errors='coerce').dropna()
                        if len(_vals):
                            _obs_median_native[f'te{_i+1}'] = float(_vals.median())
        except Exception as _e:
            log(f'Semi-quant: could not read proxy medians ({_e})')

        def _aq(row, idx):
            """aq_conc in ppb for driprates().
            In semi mode: back-calculate aq_conc (ppb) from the observed proxy median
            using the correct Kp-based kinetic model:
              h(V) = Kp*(XF+XL)*(aq/ca)*400000*K_e * (1 - nS*exp(-Kd*tau))
            Solve for aq_ppb so that h(V_ref) = obs_ppm at a reference drip rate.
            """
            if _analysis_mode == 'semi':
                obs_native = _obs_median_native.get(row)
                if obs_native and obs_native > 0:
                    kp_  = float(params.get(row + '_Kp', 1))
                    if kp_ < 0:
                        kp_ = {'Co': 4.4, 'Ni': 1.1, 'Cu': 44}.get(
                            params.get(row + '_elem', ''), 1)
                    kd_  = float(np.exp(float(params[row + '_Kd_mn'])))
                    ke_  = float(params.get(row + '_K_e', 1))
                    xf_  = float(params.get(row + '_F', 0.01))
                    xl_  = max(float(params[row + '_labile']) if params.get(row + '_labile')
                               else 1.0 - float(params.get(row + '_InertF', 0))
                                        - float(params.get(row + '_F', 0)), 0.01)
                    nS_  = 1.0 - xf_
                    # Reference drip rate for the kinetic inversion (10 drips/min default)
                    tau_ = 6.0   # seconds (= 60/10)
                    E1_  = float(np.exp(-kd_ * tau_))
                    kin_factor = 1.0 - nS_ * E1_   # fraction of K0 reached at V_ref

                    # Convert obs to ppm
                    _te_entries = params.get('te_list') or []
                    _unt = _te_entries[idx].get('unit', 'ppm') \
                           if idx < len(_te_entries) else params.get(row + '_unit', 'ppm')
                    obs_ppm  = obs_native * _unit_to_ppb(_unt) / 1000.0
                    ca_ppb_  = float(params['ca_conc']) * _unit_to_ppb(params.get('ca_unit','ppb'))

                    # Solve: obs = Kp*(XF+XL)*(aq/ca)*400000*K_e * kin_factor
                    # => aq = obs * ca / (Kp*(XF+XL)*400000*K_e*kin_factor)
                    denom = max(kp_ * (xf_ + xl_) * CA_CALCITE_PPM * ke_ * kin_factor, 1e-20)
                    implied  = obs_ppm * ca_ppb_ / denom
                    log(f'Semi-quant {row}: obs={obs_native:.4f} (native), '
                        f'obs_ppm={obs_ppm:.4f}, ca_ppb={ca_ppb_:.1f}, '
                        f'Kp={kp_}, kin_factor={kin_factor:.4f}, '
                        f'implied aq_conc={implied:.4f} ppb')
                    return implied
                log(f'Semi-quant {row}: no proxy median — using entered aq_conc')
                return float(params[row + '_aq_conc']) * _unit_to_ppb(params.get(row + '_aq_unit', 'ppb'))
            return float(params[row + '_aq_conc']) * _unit_to_ppb(params.get(row + '_aq_unit', 'ppb'))

        def make_te(pdist, row, idx):
            _kp_val = float(params[row + '_Kp'])
            _elem   = params[row + '_elem']
            # Resolve Kp=-1 (theoretical): model.py kp_theory only supports Cu, Ni, Co
            _KP_THEORY_ELEMS = {'Cu', 'Ni', 'Co'}
            if _kp_val == -1:
                if _elem not in _KP_THEORY_ELEMS:
                    # No theoretical Kp — use THEO_KP empirical values or raise error
                    _empirical = {'Co': 4.4, 'Ni': 1.1, 'Cu': 44}
                    if _elem in _empirical:
                        _kp_val = _empirical[_elem]
                        log(f'  Kp=-1 for {_elem}: resolved to empirical {_kp_val}')
                    else:
                        raise ValueError(
                            f'Kp=-1 (theoretical) is not supported for {_elem}. '
                            f'Supported elements: {", ".join(sorted(_KP_THEORY_ELEMS))}. '
                            f'Please enter a numeric Kp value for {_elem}.')

            te = {
                'elem':    _elem,
                'mol_wt':  float(params[row + '_mol_wt']),
                'Kp':      _kp_val,
                'Temp_C':  float(params['temp_C']),
                'F':       float(params[row + '_F']),
                'InertF':  float(params[row + '_InertF']),
                'aq_conc': np.float64(_aq(row, idx)),
                'ca_conc': np.float64(
                    float(params['ca_conc']) * _unit_to_ppb(params.get('ca_unit', 'ppb'))
                ),
                'PDist':   pdist,
            }

            # Attach concentration priors if supplied (stochastic mode)
            try:
                from concentration_prior import ConcentrationPrior
                if 'ca_prior_mu_ln' in params and 'ca_prior_sigma_ln' in params:
                    te['ca_prior'] = ConcentrationPrior(
                        mu_ln    = float(params['ca_prior_mu_ln']),
                        sigma_ln = float(params['ca_prior_sigma_ln']),
                        unit     = 'ppb',
                        source   = params.get('ca_prior_source', 'manual'),
                    )
                _aq_mu_key  = f'{row}_prior_mu_ln'
                _aq_sig_key = f'{row}_prior_sigma_ln'
                if _aq_mu_key in params and _aq_sig_key in params:
                    te['aq_prior'] = ConcentrationPrior(
                        mu_ln    = float(params[_aq_mu_key]),
                        sigma_ln = float(params[_aq_sig_key]),
                        unit     = 'ppb',
                        source   = params.get(f'{row}_prior_source', 'manual'),
                    )
            except ImportError:
                pass  # concentration_prior.py not available — scalar only

            # Log stochastic status
            if 'ca_prior' in te:
                log(f'  ├─ [Ca_aq] STOCHASTIC: µ_ln={te["ca_prior"].mu_ln:.3f}, '
                    f'σ_ln={te["ca_prior"].sigma_ln:.3f}, '
                    f'mean={te["ca_prior"].mean_ppb:.1f} ppb')
            else:
                log(f'  ├─ [Ca_aq] fixed: {te["ca_conc"]:.1f} ppb')
            if 'aq_prior' in te:
                log(f'  ├─ [TE_aq] STOCHASTIC: µ_ln={te["aq_prior"].mu_ln:.3f}, '
                    f'σ_ln={te["aq_prior"].sigma_ln:.3f}, '
                    f'mean={te["aq_prior"].mean_ppb:.3f} ppb')
            else:
                log(f'  ├─ [TE_aq] fixed: {te["aq_conc"]:.4f} ppb')

            return te

        # ── Run inline diagnostic test before the real driprates call ─────
        try:
            _test_te = make_te(PDist_TEs[0], 'te1', 0)
            _inline_test(PDist_TEs[0], _test_te,
                         float(params['te1_Kd_mn']), float(params['te1_Kd_sd']),
                         float(params.get('te1_K_e', 1)))
        except Exception as _e:
            log(f'  [TEST] Exception: {_e}')
            log(traceback.format_exc())

        # Loop over all TE PDists
        V_pdf_list = []
        for _i, _pd in enumerate(PDist_TEs):
            _rk = f'te{_i+1}'
            _pct = 65 + int(15 * _i / len(PDist_TEs))
            _set_stage(f'Computing drip rate PDF ({_rk.upper()})', _pct)
            log(f'Computing drip rate PDF for {_rk.upper()} …')
            _TE = make_te(_pd, _rk, _i)
            _Ke = float(params.get(f'{_rk}_K_e', params.get(f'K_e{_i+1}', 1)))

            # ── Diagnostic: log every parameter entering driprates() ──────
            _diag_kd_mn = float(params[f'{_rk}_Kd_mn'])
            _diag_kd_sd = float(params[f'{_rk}_Kd_sd'])
            _diag_kd    = float(np.exp(_diag_kd_mn))
            _diag_xl    = 1.0 - _TE['F'] - _TE['InertF']
            _diag_xf    = _TE['F']
            _diag_aq    = _TE['aq_conc']  # ppb
            _diag_ca    = _TE['ca_conc']  # ppb
            _diag_kp    = _TE['Kp']
            _diag_nS    = 1.0 - _diag_xf
            _diag_K0    = _diag_kp * (_diag_xf + _diag_xl) * (_diag_aq / _diag_ca) * CA_CALCITE_PPM * _Ke
            _diag_tau   = 6.0  # 10 drips/min reference
            _diag_E1    = float(np.exp(-_diag_kd * _diag_tau))
            _diag_kin   = _diag_K0 * (1 - _diag_nS * _diag_E1)
            log(f'  ├─ elem={_TE["elem"]}  mol_wt={_TE["mol_wt"]}  Kp={_diag_kp}  T={_TE["Temp_C"]}°C')
            log(f'  ├─ Kd_mn={_diag_kd_mn:.4f}  Kd_sd={_diag_kd_sd:.4f}  Kd={_diag_kd:.6f}  K_e={_Ke}')
            log(f'  ├─ XF={_diag_xf}  XI={_TE["InertF"]}  XL={_diag_xl:.4f}  nS={_diag_nS}')
            log(f'  ├─ aq_conc={_diag_aq:.4f} ppb  ca_conc={_diag_ca:.1f} ppb')
            log(f'  ├─ aq/ca ratio = {_diag_aq/_diag_ca:.6e}')
            log(f'  ├─ K0_equilibrium = {_diag_K0:.6f} ppm (using Kp={_diag_kp})')
            log(f'  ├─ pred_kinetic @10/min = {_diag_kin:.6f} ppm (E1={_diag_E1:.6f})')
            # PDist diagnostics — dump actual attributes since BayProX naming is unknown
            try:
                _pd_attrs = [a for a in dir(_pd) if not a.startswith('_')]
                log(f'  ├─ PDist type={type(_pd).__name__}, attrs={_pd_attrs}')
                for _attr in _pd_attrs:
                    try:
                        _val = getattr(_pd, _attr)
                        if callable(_val):
                            continue
                        if isinstance(_val, np.ndarray):
                            log(f'  ├─ PDist.{_attr}: ndarray shape={_val.shape}, '
                                f'dtype={_val.dtype}, range=[{np.nanmin(_val):.6g}, {np.nanmax(_val):.6g}], '
                                f'finite={np.count_nonzero(np.isfinite(_val))}/{_val.size}')
                        elif isinstance(_val, (int, float, str, bool)):
                            log(f'  ├─ PDist.{_attr} = {_val}')
                        elif isinstance(_val, (list, tuple)) and len(_val) < 10:
                            log(f'  ├─ PDist.{_attr} = {_val}')
                    except Exception:
                        pass
            except Exception as _de:
                log(f'  ├─ PDist inspection failed: {_de}')
            log(f'  └─ Calling driprates(Kd_mn={_diag_kd_mn}, Kd_sd={_diag_kd_sd}, K_e={_Ke})')

            _vpdf, V_age, V_span = driprates(
                _diag_kd_mn, _diag_kd_sd,
                _Ke, TE=_TE, calib=False)

            # ── Diagnostic: inspect output ────────────────────────────────
            _vpdf_nz  = np.count_nonzero(_vpdf)
            _vpdf_fin = np.count_nonzero(np.isfinite(_vpdf))
            _vpdf_nan = np.count_nonzero(np.isnan(_vpdf))
            log(f'  driprates() returned: V_pdf shape={_vpdf.shape}, '
                f'nonzero={_vpdf_nz}/{_vpdf.size}, '
                f'NaN={_vpdf_nan}, finite={_vpdf_fin}')
            log(f'  V_span: [{V_span[0]:.4f}, {V_span[-1]:.4f}] ({len(V_span)} pts)')
            log(f'  V_age:  [{V_age[0]:.1f}, {V_age[-1]:.1f}] ({len(V_age)} pts)')
            # Per-column diagnostics: how many timesteps have any finite signal?
            _col_sums = np.nansum(_vpdf, axis=0)
            _col_finite = np.count_nonzero(np.isfinite(_col_sums) & (_col_sums > 0))
            log(f'  Timesteps with finite nonzero PDF: {_col_finite}/{_vpdf.shape[1]}')
            if _col_finite == 0:
                # Extra debug: check if the model can predict anything in the PDist range
                log(f'  ⚠ ALL TIMESTEPS ZERO/NaN — diagnosing forward model vs PDist mismatch:')
                # K_0 uses Kp; kinetic range is [K_0*XF, K_0]
                log(f'    K0 = {_diag_K0:.6f} ppm (Kp={_diag_kp}, equilibrium)')
                log(f'    Model h(v) range: [{_diag_K0*_diag_xf:.6f}, {_diag_K0:.6f}] ppm')
                log(f'    (fast-only limit to full equilibrium)')

            V_pdf_list.append(np.nan_to_num(_vpdf, nan=0.0))
        V_pdf_TE1 = V_pdf_list[0]
        V_pdf_TE2 = V_pdf_list[1] if len(V_pdf_list) > 1 else None

        # ── Kd sensitivity: run at Kd_mn ± Kd_sd (near-deterministic) ──
        _kd_sens_lo_pdfs = []
        _kd_sens_hi_pdfs = []
        for _i, _pd in enumerate(PDist_TEs):
            _rk = f'te{_i+1}'
            _TE_s = make_te(_pd, _rk, _i)
            _Ke   = float(params.get(f'{_rk}_K_e', params.get(f'K_e{_i+1}', 1)))
            _kmn  = float(params[f'{_rk}_Kd_mn'])
            _ksd  = float(params[f'{_rk}_Kd_sd'])
            _eps  = 0.001  # near-zero sd → near-deterministic
            _vlo, _, _ = driprates(_kmn - _ksd, _eps, _Ke, TE=_TE_s, calib=False)
            _vhi, _, _ = driprates(_kmn + _ksd, _eps, _Ke, TE=_TE_s, calib=False)
            _kd_sens_lo_pdfs.append(np.nan_to_num(_vlo, nan=0.0))
            _kd_sens_hi_pdfs.append(np.nan_to_num(_vhi, nan=0.0))
        log('Kd sensitivity passes complete')

        # ── Kd sensitivity joint PDFs ───────────────────────────────────────
        def _joint_and_median(pdf_list, rsw):
            n = len(pdf_list)
            if n == 1:
                jp = pdf_list[0].copy()
                jp[~np.isfinite(jp) | (jp < 0)] = 0.
            else:
                jp = pdf_list[0].copy()
                for _v in pdf_list[1:]: jp = jp * _v
                jp[~np.isfinite(jp) | (jp < 0)] = 0.
                jp = np.power(jp, 1.0 / n)
            jp = np.nan_to_num(jp, nan=0.0)
            C = (jp.T * rsw).sum(axis=1)
            C[(C == 0) | ~np.isfinite(C)] = 1.0
            jp = jp / C
            med = np.zeros(jp.shape[1])
            for _i in range(jp.shape[1]):
                cdf = np.cumsum(rsw * jp[:, _i])
                if cdf[-1] <= 0 or not np.isfinite(cdf[-1]):
                    continue  # leave med[_i] = 0
                cdf /= cdf[-1]
                med[_i] = float(interp1d(cdf, V_span, kind='linear',
                    bounds_error=False, fill_value=(V_span[0], V_span[-1]))(0.5))
            return med

        # ── 7. Joint PDF ────────────────────────────────────────────────────
        _set_stage('Building joint PDF', 82)
        V_rsw = 0.5 * np.r_[
            V_span[1]  - V_span[0],
            V_span[2:] - V_span[:-2],
            V_span[-1] - V_span[-2]
        ]
        # Sensitivity medians need V_rsw — compute now
        _kd_med_lo = _joint_and_median(_kd_sens_lo_pdfs, V_rsw)
        _kd_med_hi = _joint_and_median(_kd_sens_hi_pdfs, V_rsw)

        _n_te_pdf = len(V_pdf_list)
        if _n_te_pdf == 1:
            log('Joint PDF uses single TE')
            V_pdf = V_pdf_list[0].copy()
            V_pdf[~np.isfinite(V_pdf) | (V_pdf < 0.)] = 0.
        else:
            log(f'Joint PDF: geometric mean of {_n_te_pdf} TEs')
            _prod = V_pdf_list[0].copy()
            for _v in V_pdf_list[1:]:
                _prod = _prod * _v
            _prod[~np.isfinite(_prod) | (_prod < 0.)] = 0.
            V_pdf = np.power(_prod, 1.0 / _n_te_pdf)
        V_pdf = np.nan_to_num(V_pdf, nan=0.0)
        C = (V_pdf.T * V_rsw).sum(axis=1)
        # Guard against zero-sum or NaN columns (PDF collapsed at some timesteps)
        C[(C == 0) | ~np.isfinite(C)] = 1.0
        V_pdf = V_pdf / C

        # ── 8. Summary statistics ───────────────────────────────────────────
        _set_stage('Extracting percentiles', 88)
        pcs = [5, 10, 25, 50, 75, 90, 95]
        V_pcs = {p: np.zeros(len(V_age)) for p in pcs}
        _n_zero_cols = 0
        for i in range(V_pdf.shape[1]):
            cdf = np.cumsum(V_rsw * V_pdf[:, i])
            if cdf[-1] <= 0 or not np.isfinite(cdf[-1]):
                _n_zero_cols += 1
                continue  # leave V_pcs[p][i] = 0; flagged below
            cdf /= cdf[-1]
            f_ = interp1d(cdf, V_span, kind='linear',
                          bounds_error=False, fill_value=(V_span[0], V_span[-1]))
            for p in pcs:
                V_pcs[p][i] = f_(p / 100.)
        if _n_zero_cols > 0:
            log(f'Warning: {_n_zero_cols}/{V_pdf.shape[1]} timesteps had zero PDF '
                f'— drip rate could not be estimated. Check Kd, aq_conc and data units.')

        # ── Semi-quant normalisation ─────────────────────────────────────────
        if _analysis_mode == 'semi':
            _set_stage('Normalising (semi-quantitative)', 90)
            _anchor = params.get('semi_anchor', '')
            _ref_min = params.get('semi_ref_min', '')
            _ref_max = params.get('semi_ref_max', '')

            # Determine normalisation scalar
            if _anchor and str(_anchor).strip():
                # User supplied an anchor drip rate → use as scale
                _scale = float(_anchor)
                log(f'Semi-quant: anchoring to {_scale} drips/min')
            elif _ref_min and _ref_max and str(_ref_min).strip() and str(_ref_max).strip():
                # Reference period: mean of pc50 within that age window
                _rlo, _rhi = float(_ref_min), float(_ref_max)
                _mask = (V_age >= min(_rlo,_rhi)) & (V_age <= max(_rlo,_rhi))
                if _mask.sum() > 0:
                    _scale = float(np.nanmean(V_pcs[50][_mask]))
                    log(f'Semi-quant: reference period mean = {_scale:.4f} → 100%')
                else:
                    _scale = float(np.nanmax(V_pcs[95]))
                    log('Semi-quant: reference period empty, falling back to record max')
            else:
                # Default: record maximum of pc95 = 100%
                _scale = float(np.nanmax(V_pcs[95]))
                log(f'Semi-quant: normalising to record max ({_scale:.4f})')

            if _scale > 0:
                for _p in V_pcs:
                    V_pcs[_p] = V_pcs[_p] / _scale * 100.0
                _kd_med_lo = _kd_med_lo / _scale * 100.0
                _kd_med_hi = _kd_med_hi / _scale * 100.0
            log('Semi-quant normalisation applied.')
        else:
            _scale = None  # full mode — no normalisation

        # ── 9. Realisations ─────────────────────────────────────────────────
        _gen_real = params.get('generate_realisations', True)
        if str(_gen_real).lower() not in ('false', '0', 'no'):
            n_real = int(params.get('n_realisations', 1000))
            _set_stage(f'Sampling {n_real} realisations', 92)
            log(f'Sampling {n_real} realisations for RQA …')
            _seed = int(params.get('rng_seed', 42))
            rng = np.random.default_rng(seed=_seed)
            realisations = np.zeros((n_real, len(V_age)))
            for i in range(V_pdf.shape[1]):
                cdf = np.cumsum(V_rsw * V_pdf[:, i])
                cdf /= cdf[-1]
                u = rng.uniform(0., 1., n_real)
                f_ = interp1d(cdf, V_span, kind='linear',
                              bounds_error=False,
                              fill_value=(V_span[0], V_span[-1]))
                realisations[:, i] = f_(u)
            # Semi-quant: normalise each realisation independently
            if _analysis_mode == 'semi' and _scale and _scale > 0:
                realisations = realisations / _scale * 100.0
                log('Semi-quant: ensemble normalisation applied to realisations')
        else:
            realisations = None
            log('Realisations skipped (generate_realisations=False)')

        # ── 10. Apply hiatus exclusion zones ────────────────────────────────
        _hiatus_zones = params.get('hiatus_zones', [])
        if _hiatus_zones:
            _n_masked = 0
            for _hz in _hiatus_zones:
                _hf = float(_hz.get('from', 0))
                _ht = float(_hz.get('to', 0))
                _hlo, _hhi = min(_hf, _ht), max(_hf, _ht)
                _hmask = (V_age >= _hlo) & (V_age <= _hhi)
                _n_this = int(_hmask.sum())
                _n_masked += _n_this
                for _p in V_pcs:
                    V_pcs[_p][_hmask] = np.nan
                _kd_med_lo[_hmask] = np.nan
                _kd_med_hi[_hmask] = np.nan
                V_pdf[:, _hmask] = 0.0  # zero out PDF in hiatus zones
                if realisations is not None:
                    realisations[:, _hmask] = np.nan
                log(f'Hiatus zone {_hlo:.0f}–{_hhi:.0f} yrs BP: masked {_n_this} timesteps')
            log(f'Total hiatus-masked timesteps: {_n_masked}')

        # ── 11. Save outputs ────────────────────────────────────────────────
        _set_stage('Saving outputs', 96)

        # Summary CSV
        summary_path = os.path.join(OUTPUT_FOLDER, 'drip_rate_summary.csv')
        df_out = pd.DataFrame({'age': V_age})
        for p in pcs:
            df_out[f'pc{p:02d}'] = V_pcs[p]
        df_out.to_csv(summary_path, index=False)
        log('Saved drip_rate_summary.csv')

        # Realisations CSV (only if generated)
        if realisations is not None:
            real_path = os.path.join(OUTPUT_FOLDER, 'drip_rate_realisations.csv')
            header = 'age,' + ','.join([f'r{j}' for j in range(n_real)])
            out_arr = np.vstack([V_age, realisations]).T
            np.savetxt(real_path, out_arr, delimiter=',',
                       header=header, comments='')
            log('Saved drip_rate_realisations.csv')

        # Plot data JSON (for browser chart)
        chart_path = os.path.join(OUTPUT_FOLDER, 'chart_data.json')
        # Smart & Friedrich classification stats (full mode only)
        # Mean and CV computed from median time series and IQR proxy for spread
        _med   = np.array(V_pcs[50])
        _p05   = np.array(V_pcs[5])
        _p95   = np.array(V_pcs[95])
        _p25   = np.array(V_pcs[25])
        _p75   = np.array(V_pcs[75])
        _mu    = float(np.nanmean(_med))
        _sd    = float(np.nanstd(_med))
        # In semi mode percentiles are in %, so CV is still valid but
        # absolute mean is not interpretable — omit sf block
        _cv    = float(_sd / _mu) if _mu > 0 and _analysis_mode == 'full' else 0.0
        # Per-timestep CV envelope using pc05/pc95 spread / pc50
        _cv_lo = float(np.nanmean((_med - _p05) / np.where(_med>0, _med, np.nan)))
        _cv_hi = float(np.nanmean((_p95 - _med) / np.where(_med>0, _med, np.nan)))
        # Temporal quantiles of pc50 for spread on scatter
        _mu_lo = float(np.nanpercentile(_med, 25))
        _mu_hi = float(np.nanpercentile(_med, 75))
        # NaN → None for JSON serialisation
        def _to_json_list(arr):
            return [None if (isinstance(v, float) and not np.isfinite(v)) else float(v)
                    for v in arr]

        chart_data = {
            'mode':  _analysis_mode,
            'age':   V_age.tolist(),
            'pc05':  _to_json_list(V_pcs[5]),
            'pc25':  _to_json_list(V_pcs[25]),
            'pc50':  _to_json_list(V_pcs[50]),
            'pc75':  _to_json_list(V_pcs[75]),
            'pc95':  _to_json_list(V_pcs[95]),
            'kd_lo': _to_json_list(_kd_med_lo),
            'kd_hi': _to_json_list(_kd_med_hi),
            'hiatus_zones': [{'from': float(z.get('from',0)), 'to': float(z.get('to',0))}
                             for z in _hiatus_zones] if _hiatus_zones else [],
            'sf': None if _analysis_mode == 'semi' else {
                'mean': _mu, 'cv': _cv,
                'mean_lo': _mu_lo, 'mean_hi': _mu_hi,
                'cv_lo': max(0, _cv - _cv_lo),
                'cv_hi': _cv + _cv_hi,
            } if _analysis_mode == 'full' else None,
        }
        with open(chart_path, 'w') as f:
            json.dump(chart_data, f)

        _outputs = ['drip_rate_summary.csv']
        if realisations is not None:
            _outputs.append('drip_rate_realisations.csv')

        # ── Run ID (timestamp-based) ────────────────────────────────
        # ── Run ID (timestamp-based) ────────────────────────────────
        _run_id = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        _run_id_path = os.path.join(OUTPUT_FOLDER, 'run_id.txt')
        with open(_run_id_path, 'w') as f:
            f.write(_run_id)
        _outputs.append('run_id.txt')
        log(f'Run ID: {_run_id}')

        # ── Age model JSON (for browser chart) ──────────────────────
        try:
            _age_data = {'depth': [], 'age_median': [], 'age_lo': None, 'age_hi': None,
                         'dated_depth': [], 'dated_age': [], 'dated_err': []}
            _pd0 = PDist_TEs[0]
            _age_data['age_median'] = _pd0.calbp.tolist()

            # ── Build depth grid ─────────────────────────────────────
            # Strategy: read depth_age.csv dated points, interpolate depth=f(age)
            _da_csv = os.path.join(UPLOAD_FOLDER, 'depth_age.csv')
            _da_depths, _da_ages, _da_errs = [], [], []
            if os.path.isfile(_da_csv):
                _dadf = pd.read_csv(_da_csv)
                _dcol_da = params.get('col_depth', '')
                _acol_da = params.get('col_age', '')
                _ecol_da = params.get('col_age_err', '')
                log(f'Age model: reading depth_age.csv, col_depth={_dcol_da!r}, col_age={_acol_da!r}, col_err={_ecol_da!r}')
                log(f'Age model: available columns: {list(_dadf.columns)}')

                if _dcol_da and _dcol_da in _dadf.columns and _acol_da and _acol_da in _dadf.columns:
                    _da_d = pd.to_numeric(_dadf[_dcol_da], errors='coerce')
                    _da_a = pd.to_numeric(_dadf[_acol_da], errors='coerce')
                    _valid = _da_d.notna() & _da_a.notna()
                    _da_depths = _da_d[_valid].tolist()
                    _da_ages   = _da_a[_valid].tolist()
                    if _ecol_da and _ecol_da in _dadf.columns:
                        _da_errs = pd.to_numeric(_dadf[_ecol_da], errors='coerce')[_valid].tolist()

                    # Interpolate depth = f(age) using the dated points
                    if len(_da_depths) >= 2:
                        from scipy.interpolate import interp1d as _interp
                        _sort = np.argsort(_da_ages)
                        _f_depth = _interp(
                            np.array(_da_ages)[_sort], np.array(_da_depths)[_sort],
                            kind='linear', bounds_error=False,
                            fill_value=(_da_depths[_sort[0]], _da_depths[_sort[-1]]))
                        _age_data['depth'] = _f_depth(_pd0.calbp).tolist()
                        log(f'Age model: interpolated depth grid from {len(_da_depths)} dated points')

            # Fallback: try TE CSV depth column
            if not _age_data['depth']:
                _depth_csv = os.path.join(UPLOAD_FOLDER, 'trace_elem1.csv')
                if os.path.isfile(_depth_csv):
                    _tdf = pd.read_csv(_depth_csv)
                    # Try multiple strategies to find depth column
                    _dcol = ''
                    _te_entries = params.get('te_list') or []
                    if _te_entries and isinstance(_te_entries[0], dict):
                        _dcol = _te_entries[0].get('col_depth', '')
                    if not _dcol:
                        # Guess from column names
                        for _guess in ['depth','Depth','DEPTH','dist','Dist','distance']:
                            if _guess in _tdf.columns:
                                _dcol = _guess; break
                    if _dcol and _dcol in _tdf.columns:
                        _depths = pd.to_numeric(_tdf[_dcol], errors='coerce').dropna().values
                        if len(_depths) > 0:
                            _dmin, _dmax = float(_depths.min()), float(_depths.max())
                            _age_data['depth'] = np.linspace(_dmin, _dmax, len(_pd0.calbp)).tolist()
                            log(f'Age model: linear depth grid from TE CSV [{_dmin:.1f}, {_dmax:.1f}]')

            # Last fallback: sequential indices
            if not _age_data['depth']:
                _age_data['depth'] = list(range(len(_pd0.calbp)))
                log('Age model: using sequential indices for depth (no depth column found)')

            # Dated points for overlay
            _age_data['dated_depth'] = _da_depths
            _age_data['dated_age']   = _da_ages
            _age_data['dated_err']   = _da_errs

            _am_path = os.path.join(OUTPUT_FOLDER, 'age_model.json')
            with open(_am_path, 'w') as f:
                json.dump(_age_data, f)
            _outputs.append('age_model.json')
            log(f'Saved age_model.json ({len(_age_data["depth"])} depth pts, {len(_da_depths)} dated pts)')
            # Also save as CSV
            _am_csv = os.path.join(OUTPUT_FOLDER, 'age_model.csv')
            pd.DataFrame({'depth': _age_data['depth'],
                          'age_yBP': _age_data['age_median']}).to_csv(_am_csv, index=False)
            _outputs.append('age_model.csv')
            log('Saved age_model.csv')
        except Exception as _e:
            log(f'Could not save age model: {_e}')
            log(traceback.format_exc())

        # ── PDF heatmap JSON (subsampled for browser) ────────────────
        try:
            # Auto-crop V range to where density exists (±5% padding)
            _row_sums = V_pdf.sum(axis=1)
            _total_density = _row_sums.sum()
            if _total_density > 0:
                _cumsum = np.cumsum(_row_sums) / _total_density
                _v_lo = max(0, np.searchsorted(_cumsum, 0.001) - 2)
                _v_hi = min(len(_cumsum) - 1, np.searchsorted(_cumsum, 0.999) + 2)
                # Ensure at least 20 bins
                if _v_hi - _v_lo < 20:
                    _mid = (_v_lo + _v_hi) // 2
                    _v_lo = max(0, _mid - 10)
                    _v_hi = min(len(_cumsum) - 1, _mid + 10)
                log(f'PDF heatmap: V density in bins {_v_lo}–{_v_hi} '
                    f'(V={V_span[_v_lo]:.2f}–{V_span[_v_hi]:.2f} drips/min)')
            else:
                _v_lo, _v_hi = 0, V_pdf.shape[0] - 1

            _cropped_pdf = V_pdf[_v_lo:_v_hi+1, :]
            _cropped_v   = V_span[_v_lo:_v_hi+1]

            # Subsample for manageable browser rendering
            _max_v = 200   # max V bins
            _max_t = 500   # max time steps
            _v_step = max(1, _cropped_pdf.shape[0] // _max_v)
            _t_step = max(1, _cropped_pdf.shape[1] // _max_t)
            _sub_pdf = _cropped_pdf[::_v_step, ::_t_step]
            _sub_v   = _cropped_v[::_v_step]
            _sub_age = V_age[::_t_step]
            _hm_data = {
                'V_pdf':  _sub_pdf.tolist(),
                'V_span': _sub_v.tolist(),
                'ages':   _sub_age.tolist(),
            }
            _hm_path = os.path.join(OUTPUT_FOLDER, 'pdf_heatmap.json')
            with open(_hm_path, 'w') as f:
                json.dump(_hm_data, f)
            _outputs.append('pdf_heatmap.json')
            log(f'Saved pdf_heatmap.json ({_sub_pdf.shape[0]}×{_sub_pdf.shape[1]} subsampled, '
                f'V range: {_sub_v[0]:.2f}–{_sub_v[-1]:.2f})')
        except Exception as _e:
            log(f'Could not save PDF heatmap: {_e}')

        # ── Input summary CSV ────────────────────────────────────────
        try:
            _summary_rows = [
                ('run_id', _run_id),
                ('timestamp', datetime.datetime.now().isoformat()),
                ('analysis_mode', _analysis_mode),
                ('station_name', params.get('station_name', '')),
                ('cave_temperature_C', params.get('temp_C', '')),
                ('global_drip_rate', params.get('global_drip_rate', '')),
                ('da_depth_unit', params.get('da_depth_unit', 'cm')),
                ('te_depth_unit', params.get('te_depth_unit', 'cm')),
                ('ca_conc', params.get('ca_conc', '')),
                ('ca_unit', params.get('ca_unit', '')),
                ('calage_min', params.get('calage_min', '')),
                ('calage_max', params.get('calage_max', '')),
                ('n_realisations', params.get('n_realisations', '')),
                ('rng_seed', params.get('rng_seed', '')),
                ('outlier_win_size', params.get('outlier_win_size', '')),
                ('v_max', params.get('v_max', 100)),
                ('v_res', params.get('v_res', 5000)),
                ('hiatus_zones', json.dumps(params.get('hiatus_zones', []))),
            ]
            # Per-TE params
            for _i in range(len(params.get('te_list', [{}]))):
                _rk = f'te{_i+1}'
                for _k in ['elem','mol_wt','Kp','Kd_mn','Kd_sd','K_e','cal_pct',
                            'F','InertF','labile','aq_conc','aq_unit']:
                    _summary_rows.append((f'{_rk}_{_k}', params.get(f'{_rk}_{_k}', '')))
                # Stochastic priors if present
                for _pk in ['prior_mu_ln', 'prior_sigma_ln', 'prior_source']:
                    _key = f'{_rk}_{_pk}'
                    if _key in params:
                        _summary_rows.append((_key, params[_key]))
            for _pk in ['ca_prior_mu_ln', 'ca_prior_sigma_ln', 'ca_prior_source']:
                if _pk in params:
                    _summary_rows.append((_pk, params[_pk]))

            _is_path = os.path.join(OUTPUT_FOLDER, 'input_summary.csv')
            pd.DataFrame(_summary_rows, columns=['parameter', 'value']).to_csv(_is_path, index=False)
            _outputs.append('input_summary.csv')
            log('Saved input_summary.csv')
        except Exception as _e:
            log(f'Could not save input summary: {_e}')

        # Include ProxyRecord.pkl in downloads if present
        _pkl_path = os.path.join(OUTPUT_FOLDER, 'ProxyRecord.pkl')
        if os.path.isfile(_pkl_path):
            _outputs.append('ProxyRecord.pkl')

        run_state['outputs'] = _outputs
        _set_stage('Complete', 100)
        log('✓ Run complete.')

    except Exception:
        run_state['error'] = traceback.format_exc()
        log('ERROR: ' + run_state['error'])
    finally:
        run_state['running'] = False
        run_state['done']    = True


def _set_stage(stage, progress):
    run_state['stage']    = stage
    run_state['progress'] = progress
    log(f'[{progress}%] {stage}')


@app.route('/chart_data')
def chart_data():
    path = os.path.join(OUTPUT_FOLDER, 'chart_data.json')
    if not os.path.exists(path):
        return jsonify({'error': 'No chart data yet'}), 404
    with open(path) as f:
        return jsonify(json.load(f))


# ── Entry point ──────────────────────────────────────────────────────────────
if __name__ == '__main__':
    print('\n  Drip Rate Estimator')
    print('  Open your browser at:  http://localhost:5000\n')
    app.run(debug=False, host='0.0.0.0', port=5000)
