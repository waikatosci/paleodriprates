"""
Drip Rate Estimator - Flask Web Application
============================================
Provides a browser-based interface for the speleothem drip rate
reconstruction pipeline (replaces the Excel-based workflow).
"""

import os, sys, json, threading, time, traceback, pickle, copy
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
<title>Drip Rate Estimator</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=DM+Mono:wght@400;500&family=Fraunces:ital,wght@0,300;0,600;1,300&display=swap" rel="stylesheet">
<script src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/4.4.1/chart.umd.min.js"></script>
<style>
  :root {
    --bg:       #0d1117;
    --surface:  #161b22;
    --border:   #2a3441;
    --accent:   #4cc9a0;
    --accent2:  #f7a440;
    --muted:    #6e7f8d;
    --text:     #cdd9e5;
    --texthi:   #e6edf3;
    --danger:   #e05c5c;
    --radius:   10px;
    --mono:     'DM Mono', monospace;
    --serif:    'Fraunces', Georgia, serif;
  }

  * { box-sizing: border-box; margin: 0; padding: 0; }

  body {
    background: var(--bg);
    color: var(--text);
    font-family: var(--mono);
    font-size: 13px;
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
    font-size: 22px;
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
    font-size: 11px;
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
  .field label { font-size: 11px; color: var(--muted); }
  .field input, .field select {
    background: var(--bg);
    border: 1px solid var(--border);
    border-radius: 6px;
    color: var(--texthi);
    font-family: var(--mono);
    font-size: 12px;
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
  .dropzone .dz-label { font-size: 11px; color: var(--muted); margin-bottom: 4px; }
  .dropzone .dz-name { font-size: 11px; color: var(--accent); display: none; }
  .dropzone.ready .dz-name { display: block; }
  .dropzone.ready .dz-hint { display: none; }
  .dropzone .dz-hint { font-size: 10px; color: var(--muted); }
  .dropzone input[type=file] {
    position: absolute; inset: 0; opacity: 0; cursor: pointer; width: 100%; height: 100%;
  }

  /* column selector that appears after upload */
  .col-selector { margin-top: 10px; display: none; }
  .col-selector.show { display: block; }
  .col-selector label { font-size: 10px; color: var(--muted); display: block; margin-bottom: 3px; }
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
    left: 0;
    width: 640px;
    background: #1c2635;
    border: 1px solid var(--accent);
    border-radius: 8px;
    padding: 12px 14px;
    font-size: 11px;
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
           target="_blank" title="Hartland et al., Geochemical Perspectives Letters">
          Hartland et al., Geochem. Persp. Lett.
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
        Expected format: comma-separated CSV with a header row. Depth in <strong>cm</strong>, age in <strong>years BP</strong>. Trace element concentration units are selectable below. Isotope data is <strong>optional</strong>.
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
              <label>Depth column</label>
              <select id="sel-depth_age-depth" onchange="storeCol('depth_age','depth',this.value);renderAgePlot()"></select>
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
          <div style="margin-bottom:10px">
            <label style="font-size:11px;color:var(--muted)">Depth column</label>
            <select id="te-depth-sel" style="width:100%;margin-top:4px" onchange="teDepthChanged(this.value)"></select>
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
            <div style="height:140px;margin-bottom:10px">
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

      <!-- TE1 -->
      <div class="card" id="te-param-card-1">
        <div class="card-title">⚗️ Trace Element 1 — Ni (default)</div>
        <div class="form-grid">
          <div class="field">
            <label>Element</label>
            <select id="te1_elem" onchange="updateMolWt('te1', this.value)">
              <optgroup label="d-block (well characterised)">
                <option value="Ni" selected>Ni — Nickel</option>
                <option value="Co">Co — Cobalt</option>
                <option value="Cu">Cu — Copper (high Kd, use with caution)</option>
              </optgroup>
              <optgroup label="d-block (limited characterisation)">
                <option value="V">V — Vanadium</option>
                <option value="Zn">Zn — Zinc</option>
                <option value="Cd">Cd — Cadmium</option>
              </optgroup>
              <optgroup label="p-block">
                <option value="Al">Al — Aluminium</option>
                <option value="Pb">Pb — Lead</option>
              </optgroup>
              <optgroup label="REE">
                <option value="La">La — Lanthanum</option>
                <option value="Ce">Ce — Cerium</option>
              </optgroup>
              <option value="other">Other (set mol. weight manually)</option>
            </select>
          </div>
          <div class="field">
            <label>Molecular weight (g/mol)</label>
            <input type="number" id="te1_mol_wt" value="58.693" step="0.001">
          </div>
          <div class="field has-tip">
            <div class="tip-label">
              <label>Partition coefficient Kp</label>
              <span class="tip-icon">?</span>
            </div>
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
            </div>
            <div style="display:flex; gap:6px; align-items:center">
              <input type="number" id="te1_Kp" value="1.1" step="0.01" style="flex:1"
                     oninput="updateTheoKp('te1')">
              <input type="text" id="te1_Kp_theo" readonly placeholder="theo. Kp"
                     title="Theoretical Kp from Wang & Xu (2001), calculated when Kp = −1"
                     style="flex:1; opacity:0.55; cursor:not-allowed; font-size:10px; text-align:center">
            </div>
          </div>
          <div class="field">
            <label>Mean ln(K<sub>d</sub>)</label>
            <input type="number" id="te1_Kd_mn" value="-3.908" step="0.001">
          </div>
          <div class="field">
            <label>Std dev ln(K<sub>d</sub>)</label>
            <input type="number" id="te1_Kd_sd" value="1.385" step="0.001">
          </div>
          <div class="field">
            <label>X<sub>F</sub> — Fast fraction</label>
            <input type="number" id="te1_F" value="0.01" step="0.001" min="0" max="1"
                   oninput="updateLabile('te1')">
          </div>
          <div class="field">
            <label>X<sub>I</sub> — Inert fraction</label>
            <input type="number" id="te1_InertF" value="0.1" step="0.001" min="0" max="1"
                   oninput="updateLabile('te1')">
          </div>
          <div class="field">
            <label>X<sub>L</sub> — Labile fraction (auto)</label>
            <input type="number" id="te1_labile" value="0.89" step="0.001" readonly
                   style="opacity:0.6; cursor:not-allowed">
          </div>
          <div class="field fullonly">
            <label>Aqueous concentration</label>
            <input type="number" id="te1_aq_conc" value="4.370" step="0.001">
          <select id="te1_aq_unit" style="margin-top:4px;width:100%"><option value="ppb">ppb</option><option value="ppm">ppm</option><option value="ug/g">µg/g</option><option value="mg/kg">mg/kg</option></select>
          </div>
        </div>
      </div>

      <!-- TE2 -->
      <div class="card" id="te-param-card-2">
        <div class="card-title">⚗️ Trace Element 2 — Co (default)</div>
        <div class="form-grid">
          <div class="field">
            <label>Element</label>
            <select id="te2_elem" onchange="updateMolWt('te2', this.value)">
              <optgroup label="d-block (well characterised)">
                <option value="Ni">Ni — Nickel</option>
                <option value="Co" selected>Co — Cobalt</option>
                <option value="Cu">Cu — Copper (high Kd, use with caution)</option>
              </optgroup>
              <optgroup label="d-block (limited characterisation)">
                <option value="V">V — Vanadium</option>
                <option value="Zn">Zn — Zinc</option>
                <option value="Cd">Cd — Cadmium</option>
              </optgroup>
              <optgroup label="p-block">
                <option value="Al">Al — Aluminium</option>
                <option value="Pb">Pb — Lead</option>
              </optgroup>
              <optgroup label="REE">
                <option value="La">La — Lanthanum</option>
                <option value="Ce">Ce — Cerium</option>
              </optgroup>
              <option value="other">Other (set mol. weight manually)</option>
            </select>
          </div>
          <div class="field">
            <label>Molecular weight (g/mol)</label>
            <input type="number" id="te2_mol_wt" value="58.933" step="0.001">
          </div>
          <div class="field has-tip">
            <div class="tip-label">
              <label>Partition coefficient Kp</label>
              <span class="tip-icon">?</span>
            </div>
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
            </div>
            <div style="display:flex; gap:6px; align-items:center">
              <input type="number" id="te2_Kp" value="4.4" step="0.01" style="flex:1"
                     oninput="updateTheoKp('te2')">
              <input type="text" id="te2_Kp_theo" readonly placeholder="theo. Kp"
                     title="Theoretical Kp from Wang & Xu (2001), calculated when Kp = −1"
                     style="flex:1; opacity:0.55; cursor:not-allowed; font-size:10px; text-align:center">
            </div>
          </div>
          <div class="field">
            <label>Mean ln(K<sub>d</sub>)</label>
            <input type="number" id="te2_Kd_mn" value="-5.572" step="0.001">
          </div>
          <div class="field">
            <label>Std dev ln(K<sub>d</sub>)</label>
            <input type="number" id="te2_Kd_sd" value="1.385" step="0.001">
          </div>
          <div class="field">
            <label>X<sub>F</sub> — Fast fraction</label>
            <input type="number" id="te2_F" value="0.01" step="0.001" min="0" max="1"
                   oninput="updateLabile('te2')">
          </div>
          <div class="field">
            <label>X<sub>I</sub> — Inert fraction</label>
            <input type="number" id="te2_InertF" value="0.4" step="0.001" min="0" max="1"
                   oninput="updateLabile('te2')">
          </div>
          <div class="field">
            <label>X<sub>L</sub> — Labile fraction (auto)</label>
            <input type="number" id="te2_labile" value="0.59" step="0.001" readonly
                   style="opacity:0.6; cursor:not-allowed">
          </div>
          <div class="field">
            <label>Aqueous concentration</label>
            <input type="number" id="te2_aq_conc" value="0.460" step="0.001">
          <select id="te2_aq_unit" class="fullonly" style="margin-top:4px;width:100%"><option value="ppb">ppb</option><option value="ppm">ppm</option><option value="ug/g">µg/g</option><option value="mg/kg">mg/kg</option></select>
          </div>
        </div>
      </div>

      <!-- Dynamic TE3+ param cards -->
      <div id="te-extra-param-cards"></div>

      <!-- Cave conditions -->
      <div class="card">
        <div class="card-title">🌡 Cave Conditions</div>
        <div class="form-grid">
          <div class="field">
            <label>Cave temperature (°C)</label>
            <input type="number" id="temp_C" value="16.5" step="0.1">
          </div>
          <div class="field fullonly">
            <label>Ca concentration in drip water</label>
            <input type="number" id="ca_conc" value="62000" step="100">
          <select id="ca_unit" style="margin-top:4px;width:100%"><option value="ppb">ppb</option><option value="ppm">ppm</option><option value="ug/g">µg/g</option><option value="mg/kg">mg/kg</option></select>
          </div>
        </div>
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
            <input type="checkbox" id="generate_realisations" checked
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
      </div>

      <!-- Time series panel -->
      <div id="rtab-ts">
        <div class="card">
          <div class="card-title">📈 Drip Rate vs Time</div>
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
          and thus precipitation — from stalagmite trace metals (Co, Ni). The underlying model
          links the concentration of trace elements in speleothem calcite to past drip rates
          via the kinetics of metal-organic dissociation, calibrated through monitoring data
          and refined using Monte Carlo propagation of paleo-temperature uncertainty.<br><br>
          The method was developed and validated using stalagmite HS4 from Heshang Cave, China,
          reconstructing Holocene rainfall for the Yangtze region and providing new insights
          into East Asian Summer Monsoon dynamics over the past 9,500 years. The probabilistic
          reconstruction framework — including the full ensemble of drip rate realisations
          available for download — enables recurrence quantification analysis (RQA) and other
          nonlinear time series methods.
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
            <div style="display:flex; gap:14px; align-items:flex-start; padding:10px 14px;
                        background:var(--bg); border:1px solid var(--border); border-radius:8px">
              <span style="font-size:20px">🇪🇺</span>
              <div>
                <div style="color:var(--texthi); font-weight:500">EU Horizon 2020 — Marie Skłodowska-Curie Actions</div>
                <div style="color:var(--muted); font-size:11px">Grant no. 691037 —
                  <a href="https://quest.pik-potsdam.de/" target="_blank"
                     style="color:var(--accent)">QUEST: QUantitative paleoEnvironments from SpeleoThems</a>
                </div>
              </div>
            </div>
            <div style="display:flex; gap:14px; align-items:flex-start; padding:10px 14px;
                        background:var(--bg); border:1px solid var(--border); border-radius:8px">
              <span style="font-size:20px">🇳🇿</span>
              <div>
                <div style="color:var(--texthi); font-weight:500">Te Apārangi — Royal Society of New Zealand</div>
                <div style="color:var(--muted); font-size:11px">Grant RIS-UOW1501 &nbsp;·&nbsp; Rutherford Discovery Fellowship RDF-UOW1601</div>
              </div>
            </div>
            <div style="display:flex; gap:14px; align-items:flex-start; padding:10px 14px;
                        background:var(--bg); border:1px solid var(--border); border-radius:8px">
              <span style="font-size:20px">🇳🇿</span>
              <div>
                <div style="color:var(--texthi); font-weight:500">Ministry for Business, Innovation and Employment (MBIE)</div>
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
  // Find matching nav item safely without template-literal querySelector
  document.querySelectorAll('.nav-item').forEach(n => {
    const oc = n.getAttribute('onclick') || '';
    if (oc.includes("'" + name + "'")) n.classList.add('active');
  });
  if (name === 'run') { buildSummary(); updateRuntimeEstimate(); runSanityChecks(); }
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
        populateColSelectors(key, data[key].columns);
        if (key === 'depth_age') renderAgePlot();
      }
    });
}

// ── TE row state ─────────────────────────────────────────────────────────
let teColumns = [];   // all columns from uploaded TE CSV
let teDepthCol = '';  // selected depth column
let teRowData  = [];  // [{col, unit}, ...]

const UNIT_OPTS = [
  {v:'ppb', l:'ppb (µg/L or ng/g)'},
  {v:'ppm', l:'ppm (µg/g or mg/kg)'},
  {v:'ug/g', l:'µg/g'},
  {v:'mg/kg', l:'mg/kg'},
];
const ELEM_OPTS = `
  <optgroup label="d-block (well characterised)">
    <option value="Ni">Ni — Nickel</option>
    <option value="Co">Co — Cobalt</option>
    <option value="Cu">Cu — Copper (high Kd)</option>
  </optgroup>
  <optgroup label="d-block (limited)">
    <option value="V">V — Vanadium</option>
    <option value="Zn">Zn — Zinc</option>
    <option value="Cd">Cd — Cadmium</option>
  </optgroup>
  <optgroup label="p-block">
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
  renderTERows();
  updatePreviewChart();
}

function addTERow() {
  const nonDepth = teColumns.filter(c => c !== teDepthCol);
  const used = teRowData.map(r => r.col);
  const next = nonDepth.find(c => !used.includes(c)) || nonDepth[0] || '';
  teRowData.push({col: next, unit: 'ppb'});
  renderTERows();
  renderTEParamCards();
  applyElemDetection(teRowData.length - 1);
}

function removeTERow(idx) {
  if (teRowData.length <= 1) return;
  teRowData.splice(idx, 1);
  renderTERows();
  renderTEParamCards();
}

function teSyncRow(idx) {
  const rowEl = document.querySelector(`.te-row[data-idx="${idx}"]`);
  if (!rowEl) return;
  teRowData[idx].col  = rowEl.querySelector('.te-col-sel').value;
  teRowData[idx].unit = rowEl.querySelector('.te-unit-sel').value;
  applyElemDetection(idx);
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
}

function makeTEParamCard(n) {
  const p = `te${n}`;
  return `<div class="card" id="te-param-card-${n}">
    <div class="card-title">⚗️ Trace Element ${n}</div>
    <div class="form-grid">
      <div class="field"><label>Element</label>
        <select id="${p}_elem" onchange="updateMolWt('${p}', this.value)">${ELEM_OPTS}</select></div>
      <div class="field"><label>Molecular weight (g/mol)</label>
        <input type="number" id="${p}_mol_wt" value="58.693" step="0.001"></div>
      <div class="field"><label>Partition coefficient Kp</label>
        <input type="number" id="${p}_Kp" value="-1" step="0.01" oninput="updateTheoKp('${p}')"></div>
      <div class="field"><label>Mean ln(K<sub>d</sub>)</label>
        <input type="number" id="${p}_Kd_mn" value="-3.908" step="0.001"></div>
      <div class="field"><label>Std dev ln(K<sub>d</sub>)</label>
        <input type="number" id="${p}_Kd_sd" value="1.385" step="0.001"></div>
      <div class="field"><label>X<sub>F</sub> — Fast fraction</label>
        <input type="number" id="${p}_F" value="0.01" step="0.001" min="0" max="1" oninput="updateLabile('${p}')"></div>
      <div class="field"><label>X<sub>I</sub> — Inert fraction</label>
        <input type="number" id="${p}_InertF" value="0.1" step="0.001" min="0" max="1" oninput="updateLabile('${p}')"></div>
      <div class="field"><label>X<sub>L</sub> — Labile (auto)</label>
        <input type="number" id="${p}_labile" value="0.89" step="0.001" readonly
               style="opacity:0.6;cursor:not-allowed"></div>
      <div class="field fullonly"><label>Aqueous concentration</label>
        <input type="number" id="${p}_aq_conc" value="4.370" step="0.001">
        <select id="${p}_aq_unit" style="margin-top:4px;width:100%">
          <option value="ppb">ppb</option><option value="ppm">ppm</option>
          <option value="ug/g">µg/g</option><option value="mg/kg">mg/kg</option>
        </select></div>
    </div></div>`;
}

function renderTEParamCards() {
  const n = teRowData.length;
  // TE2 built-in card: show only when ≥2 TEs
  const c2 = document.getElementById('te-param-card-2');
  if (c2) c2.style.display = n >= 2 ? '' : 'none';
  // TE3+ dynamic cards
  const container = document.getElementById('te-extra-param-cards');
  if (!container) return;
  container.innerHTML = '';
  for (let i = 2; i < n; i++) {
    container.insertAdjacentHTML('beforeend', makeTEParamCard(i + 1));
    // initialise mol wt / Kp — detect from column name, fallback to Ni
    const detectedElem = detectElemFromCol(teRowData[i]?.col) || 'Ni';
    updateMolWt(`te${i+1}`, detectedElem);
    updateLabile(`te${i+1}`);
  }
  // Re-apply detection for built-in cards (TE1, TE2) too
  for (let i = 0; i < Math.min(n, 2); i++) applyElemDetection(i);
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
    teRowData = [{col: nonDepth[0] || columns[0], unit: 'ppb'}];
    renderTERows();
    renderTEParamCards();
    // Auto-detect element for each initial row after cards are in DOM
    teRowData.forEach((_, i) => applyElemDetection(i));
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
    const guesses = {
      depth:   ['depth','Depth','DEPTH','dist','Dist'],
      age:     ['age','Age','AGE','yr_bp','yBP'],
      age_err: ['error','Error','err','Err','sigma'],
      proxy:   columns.filter(c => !c.toLowerCase().includes('depth')),
    };
    const g = (guesses[role] || []).find(h => columns.includes(h));
    if (g) sel.value = g;
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

  // Zip and filter to numeric pairs
  const pts = rawD.map((d,i) => ({d: parseFloat(d), a: parseFloat(rawA[i]),
                                   e: parseFloat(rawE[i]) || 0}))
                  .filter(p => isFinite(p.d) && isFinite(p.a));
  if (pts.length === 0) return;

  pts.sort((a,b) => a.d - b.d);
  const depths = pts.map(p=>p.d), ages = pts.map(p=>p.a);

  // Choose fit type
  const fitType = document.getElementById('age-fit-type')?.value || 'pchip';
  const fit = fitType === 'linear' ? linReg(depths, ages) : pchip(depths, ages);

  // Depth grid for smooth curve — 300 pts across extended range
  const dMin = depths[0], dMax = depths[depths.length-1];
  const dSpan = dMax - dMin || 1;
  const dExtMin = dMin - dSpan * 0.12;
  const dExtMax = dMax + dSpan * 0.12;
  const N_CURVE = 300;
  const dStep   = (dExtMax - dExtMin) / (N_CURVE - 1);
  const curvePts = Array.from({length: N_CURVE}, (_, i) => {
    const d = dExtMin + i * dStep;
    return {x: d, y: fit ? fit.predict(d) : null};
  }).filter(p => p.y !== null);

  // Split into data-range (solid) and extrapolation (dashed) segments
  const solidPts = curvePts.filter(p => p.x >= dMin && p.x <= dMax);
  const extLoPts = curvePts.filter(p => p.x <  dMin);
  const extHiPts = curvePts.filter(p => p.x >  dMax);
  // Join extrapolation segments to the nearest data-range endpoint
  if (solidPts.length) {
    if (extLoPts.length) extLoPts.push(solidPts[0]);
    if (extHiPts.length) extHiPts.unshift(solidPts[solidPts.length-1]);
  }

  // Auto-populate age range from extrapolated endpoints
  let ageMin, ageMax;
  if (fit) {
    ageMin = Math.round(Math.min(fit.predict(dExtMin), fit.predict(dExtMax)));
    ageMax = Math.round(Math.max(fit.predict(dExtMin), fit.predict(dExtMax)));
  } else {
    ageMin = Math.round(Math.min(...ages));
    ageMax = Math.round(Math.max(...ages));
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

  // Build chart datasets
  const scatterPts = pts.map(p => ({x: p.d, y: p.a}));
  const errBars    = pts.map(p => ({x: p.d, y: p.a, e: p.e}));

  const datasets = [];

  // Error bars
  if (rawE.length) {
    datasets.push({
      label: '2σ error',
      data: errBars,
      type: 'scatter',
      pointRadius: 0,
      showLine: false,
    });
  }

  // Solid fit line (within data range)
  if (solidPts.length) {
    datasets.push({
      label: fitType === 'linear' ? 'Linear fit' : 'PCHIP spline',
      data: solidPts,
      type: 'line',
      borderColor: 'rgba(76,201,160,0.7)',
      borderWidth: 1.8,
      pointRadius: 0,
      fill: false,
      tension: 0,
    });
  }

  // Dashed extrapolation below data
  if (extLoPts.length > 1) {
    datasets.push({
      label: 'Extrapolation',
      data: extLoPts,
      type: 'line',
      borderColor: 'rgba(180,180,180,0.5)',
      borderWidth: 1.5,
      borderDash: [4,3],
      pointRadius: 0,
      fill: false,
      tension: 0,
    });
  }

  // Dashed extrapolation above data
  if (extHiPts.length > 1) {
    datasets.push({
      label: 'Extrapolation (hi)',
      data: extHiPts,
      type: 'line',
      borderColor: 'rgba(180,180,180,0.5)',
      borderWidth: 1.5,
      borderDash: [4,3],
      pointRadius: 0,
      fill: false,
      tension: 0,
    });
  }

  datasets.push({
    label: 'Dated points',
    data: scatterPts,
    type: 'scatter',
    backgroundColor: 'rgba(76,201,160,0.85)',
    pointRadius: 5,
    pointHoverRadius: 7,
  });

  // Age-range shading plugin
  const shadingPlugin = {
    id: 'ageRangeShading',
    afterDraw(chart) {
      const {ctx: c, chartArea: {left,right,top,bottom}, scales: {y}} = chart;
      const mn = parseFloat(document.getElementById('calage_min')?.value);
      const mx = parseFloat(document.getElementById('calage_max')?.value);
      if (!isFinite(mn) || !isFinite(mx)) return;
      const yMn = y.getPixelForValue(mn), yMx = y.getPixelForValue(mx);
      const yTop = Math.min(yMn, yMx), yBot = Math.max(yMn, yMx);
      c.save();
      c.fillStyle = 'rgba(76,201,160,0.07)';
      c.fillRect(left, yTop, right-left, yBot-yTop);
      c.strokeStyle = 'rgba(76,201,160,0.4)';
      c.setLineDash([3,3]);
      c.lineWidth = 1;
      c.beginPath(); c.moveTo(left,yMn); c.lineTo(right,yMn); c.stroke();
      c.beginPath(); c.moveTo(left,yMx); c.lineTo(right,yMx); c.stroke();
      c.restore();
    }
  };

  agePlotChart = new Chart(ctx, {
    data: {datasets},
    options: {
      animation: false,
      responsive: true,
      maintainAspectRatio: false,
      plugins: {
        legend: {display: false},
        tooltip: {
          callbacks: {
            label: item => `depth ${item.parsed.x.toFixed(2)} cm → ${item.parsed.y.toFixed(0)} yrs BP`
          }
        }
      },
      scales: {
        x: {
          title: {display:true, text:'Depth (cm)', color:'var(--muted)', font:{size:10}},
          grid: {color:'rgba(255,255,255,0.05)'},
          ticks: {color:'var(--muted)', font:{size:10}},
        },
        y: {
          title: {display:true, text:'Age (yrs BP)', color:'var(--muted)', font:{size:10}},
          grid: {color:'rgba(255,255,255,0.05)'},
          ticks: {color:'var(--muted)', font:{size:10}},
        }
      }
    },
    plugins: [shadingPlugin],
  });
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

function updatePreviewChart() {
  // Show first configured proxy column as scatter; highlight sigma outliers
  const proxyCol = teRowData[0]?.col;
  if (!proxyCol || !teRawData[proxyCol]) return;
  const raw = teRawData[proxyCol].map(Number);
  const dep = teRawDepth;
  if (!raw.length || raw.length !== dep.length) return;

  const sigma   = parseFloat(document.getElementById('te-sigma').value) || 3;
  const targetN = parseInt(document.getElementById('te-target-n').value) || teRawN;
  const winSize = parseInt(document.getElementById('te-winsize').value) || 0;

  // Windowed (or global) sigma-clip
  const isOutlier = windowedSigmaClip(raw, sigma, winSize);

  // Block-average positions (just show midpoints)
  const step = Math.max(1, Math.round(teRawN / targetN));
  const kept = [], ktd = [], out = [], outd = [];
  raw.forEach((v, i) => {
    if (!isFinite(v)) return;
    if (isOutlier[i]) { out.push({x: dep[i], y: v}); return; }
    if (i % step === 0) ktd.push({x: dep[i], y: v});
  });

  const ctx = document.getElementById('te-preview-chart').getContext('2d');
  if (tePreviewChart) tePreviewChart.destroy();
  tePreviewChart = new Chart(ctx, {
    type: 'scatter',
    data: { datasets: [
      { label: 'Kept', data: ktd,
        backgroundColor: 'rgba(76,201,160,0.5)', pointRadius: 2, pointHoverRadius: 4 },
      { label: 'Outlier (σ-clip)', data: out,
        backgroundColor: 'rgba(255,80,80,0.7)', pointRadius: 3 },
    ]},
    options: {
      animation: false, responsive: true, maintainAspectRatio: false,
      plugins: { legend: { labels: { color:'var(--muted)', font:{size:9} } } },
      scales: {
        x: { ticks:{color:'var(--muted)',font:{size:9}}, grid:{color:'rgba(255,255,255,0.04)'},
             title:{display:true,text:'Depth (cm)',color:'var(--muted)',font:{size:9}} },
        y: { ticks:{color:'var(--muted)',font:{size:9}}, grid:{color:'rgba(255,255,255,0.04)'} }
      }
    }
  });

  const nOut = out.length;
  const nFinal = Math.round((teRawN - nOut) / Math.max(1, step));
  document.getElementById('te-row-badge').textContent =
    `${teRawN} rows → ~${nFinal} after preprocessing (${nOut} σ-outliers)`;
  updateRuntimeEstimate();
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
// BayProX dominates; drip rate step is fast (~seconds).
const BAYPROX_SEC_PER_PT_STEP = 0.000028;

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

  const bayproxSec = useCached ? 0 : nPts * nAge * BAYPROX_SEC_PER_PT_STEP * nTEs;
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
    ['Ca conc',          v('ca_conc') + ' ' + (document.getElementById('ca_unit')?.value||'ppb')],
    ['Realisations',     document.getElementById('generate_realisations')?.checked
                         ? v('n_realisations') : 'skipped'],
    ['Use cached proxy', document.getElementById('use_cached_proxy').checked ? 'Yes' : 'No'],
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
// ── Pre-run sanity checks ────────────────────────────────────────────────
const UNIT_TO_PPB = {ppb:1, ppm:1000, 'ug/g':1000, 'mg/kg':1000};

function runSanityChecks() {
  const warnings = [];
  const mode = analysisMode;

  teRowData.forEach((r, i) => {
    const n = i + 1;
    const p = `te${n}`;

    // 1. Predicted vs observed speleothem concentration
    if (mode === 'full') {
      const aqVal    = parseFloat(document.getElementById(`${p}_aq_conc`)?.value) || 0;
      const aqUnit   = document.getElementById(`${p}_aq_unit`)?.value || 'ppb';
      const aqPPB    = aqVal * (UNIT_TO_PPB[aqUnit] || 1);
      const caVal    = parseFloat(document.getElementById('ca_conc')?.value) || 0;
      const caUnit   = document.getElementById('ca_unit')?.value || 'ppb';
      const caPPB    = caVal * (UNIT_TO_PPB[caUnit] || 1);
      const kdMn     = parseFloat(document.getElementById(`${p}_Kd_mn`)?.value) || 0;
      const xl       = parseFloat(document.getElementById(`${p}_labile`)?.value) || 0;
      const dataUnit = r.unit || 'ppb';
      const dataCol  = r.col;

      if (aqPPB > 0 && caPPB > 0 && teRawData[dataCol]) {
        // Predicted speleothem conc (ppb) at median drip rate
        const kd        = Math.exp(kdMn);
        const predPPB   = kd * xl * (aqPPB / caPPB) * caPPB;  // Kd×XL×(aq/Ca)×Ca = Kd×XL×aq
        // Observed median — convert raw data to ppb
        const rawVals   = (teRawData[dataCol] || []).map(Number).filter(isFinite);
        if (rawVals.length) {
          rawVals.sort((a,b)=>a-b);
          const obsPPB  = rawVals[Math.floor(rawVals.length/2)] * (UNIT_TO_PPB[dataUnit] || 1);
          const ratio   = obsPPB / predPPB;
          if (ratio > 50 || ratio < 0.02) {
            const dir   = ratio > 1 ? 'higher' : 'lower';
            const hint  = ratio > 50
              ? 'Check aq_conc units — value may be in ppm but selector set to ppb, or ca_conc may be too high.'
              : 'Predicted concentration far exceeds data — check Kd, aq_conc, or ca_conc values.';
            warnings.push(
              `<div style="margin-bottom:10px;padding:8px 10px;background:rgba(255,160,50,0.08);
                           border-left:3px solid #ffa032;border-radius:0 4px 4px 0">
                <strong style="color:#ffa032">TE${n} (${p}_elem value) — concentration mismatch (${ratio.toFixed(0)}×)</strong><br>
                Predicted speleothem conc: <code>${predPPB.toFixed(3)} ppb</code> &nbsp;·&nbsp;
                Observed median: <code>${obsPPB.toFixed(3)} ppb</code> &nbsp;·&nbsp;
                Observed is ${ratio.toFixed(0)}× ${dir} than predicted.<br>
                <span style="color:var(--muted)">${hint}</span>
              </div>`
            );
          }
        }
      }
    }

    // 2. aq_conc unit vs data unit mismatch hint
    if (mode === 'full') {
      const aqUnit   = document.getElementById(`${p}_aq_unit`)?.value || 'ppb';
      const dataUnit = r.unit || 'ppb';
      const aqScale  = UNIT_TO_PPB[aqUnit]  || 1;
      const datScale = UNIT_TO_PPB[dataUnit] || 1;
      if (aqScale !== datScale) {
        warnings.push(
          `<div style="margin-bottom:10px;padding:8px 10px;background:rgba(100,140,255,0.08);
                       border-left:3px solid #6495ed;border-radius:0 4px 4px 0">
            <strong style="color:#6495ed">TE${n} — unit mismatch note</strong><br>
            Data column unit is <code>${dataUnit}</code> but aq_conc unit is <code>${aqUnit}</code>.
            Both are converted to ppb internally — just confirm this is intentional.
          </div>`
        );
      }
    }

    // 3. Missing aq_conc in full mode
    if (mode === 'full') {
      const aqVal = parseFloat(document.getElementById(`${p}_aq_conc`)?.value);
      if (!aqVal || aqVal <= 0) {
        warnings.push(
          `<div style="margin-bottom:10px;padding:8px 10px;background:rgba(255,80,80,0.08);
                       border-left:3px solid #ff5050;border-radius:0 4px 4px 0">
            <strong style="color:#ff5050">TE${n} — aq_conc is zero or missing</strong><br>
            Aqueous concentration must be > 0 in Full Quantification mode.
          </div>`
        );
      }
    }
  });

  // 4. Ca conc check
  if (mode === 'full') {
    const caVal  = parseFloat(document.getElementById('ca_conc')?.value) || 0;
    const caUnit = document.getElementById('ca_unit')?.value || 'ppb';
    const caMgL  = caVal * (UNIT_TO_PPB[caUnit] || 1) / 1000;  // convert ppb → mg/L
    if (caMgL > 200) {
      warnings.push(
        `<div style="margin-bottom:10px;padding:8px 10px;background:rgba(255,160,50,0.08);
                     border-left:3px solid #ffa032;border-radius:0 4px 4px 0">
          <strong style="color:#ffa032">Ca concentration unusually high (${caMgL.toFixed(0)} mg/L)</strong><br>
          Typical cave drip water Ca is 20–100 mg/L. Check units and value.
        </div>`
      );
    }
  }

  const card = document.getElementById('sanity-card');
  const body = document.getElementById('sanity-body');
  if (!card || !body) return warnings.length === 0;
  if (warnings.length === 0) {
    card.style.display = 'none';
  } else {
    body.innerHTML = warnings.join('');
    card.style.display = '';
  }
  return warnings.filter(w => w.includes('ff5050') || w.includes('ffa032')).length === 0;
                          // true = ok to run (only blue info warnings present)
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
    // TE list — one entry per configured trace element row
    te_list: teRowData.map(r => ({col_depth: teDepthCol, col_proxy: r.col, unit: r.unit})),
    // TE param cards (te1_, te2_, te3_, ... collected dynamically)
    ...(() => {
      const out = {};
      for (let i = 1; i <= Math.max(teRowData.length, 2); i++) {
        const p = `te${i}`;
        ['elem','mol_wt','Kp','Kd_mn','Kd_sd','F','InertF','aq_conc','aq_unit'].forEach(f => {
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
    ca_conc:      v('ca_conc'),
    ca_unit:      document.getElementById('ca_unit')?.value || 'ppb',
    // output
    n_realisations:       v('n_realisations'),
    rng_seed:             v('rng_seed'),
    generate_realisations: document.getElementById('generate_realisations')?.checked ?? true,
    analysis_mode:        analysisMode,
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
}

// ── Chart ────────────────────────────────────────────────────────────────────
let sfChart = null;

function switchResultTab(name) {
  ['ts','sf'].forEach(t => {
    document.getElementById(`rtab-${t}`).style.display = t===name ? '' : 'none';
    const btn = document.getElementById(`tab-${t}`);
    btn.style.borderBottomColor = t===name ? 'var(--accent)' : 'transparent';
    btn.style.color = t===name ? 'var(--texthi)' : 'var(--muted)';
  });
  if (name==='sf') loadSFChart();
}

function loadSFChart() {
  fetch('/chart_data').then(r => r.json()).then(d => {
    if (!d.sf) return;
    if (d.mode === 'semi') {
      const ctx = document.getElementById('sfChart').getContext('2d');
      if (sfChart) sfChart.destroy();
      ctx.font = '13px DM Mono, monospace';
      ctx.fillStyle = 'var(--muted)';
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
                    color:'var(--muted)', font:{size:11}},
            ticks: {color:'var(--muted)', font:{size:10}},
            grid:  {color:'rgba(255,255,255,0.05)'},
          },
          y: {
            min: 0, max: 60,
            title: {display:true, text:'Mean drip rate (drips min⁻¹)',
                    color:'var(--muted)', font:{size:11}},
            ticks: {color:'var(--muted)', font:{size:10}},
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

function loadChart() {
  fetch('/chart_data').then(r => r.json()).then(d => {
    if (d.error) return;
    const ctx = document.getElementById('dripChart').getContext('2d');
    if (chart) chart.destroy();
    chart = new Chart(ctx, {
      type: 'line',
      data: {
        labels: d.age,
        datasets: [
          {
            label: 'pc95–pc05 range',
            data: d.pc95,
            fill: '-1',
            backgroundColor: 'rgba(76,201,160,0.08)',
            borderColor: 'transparent',
            pointRadius: 0,
            tension: 0.3,
          },
          {
            label: 'pc05',
            data: d.pc05,
            borderColor: 'rgba(76,201,160,0.25)',
            borderWidth: 1,
            pointRadius: 0,
            fill: false,
            tension: 0.3,
          },
          {
            label: 'IQR (25–75)',
            data: d.pc75,
            fill: '+1',
            backgroundColor: 'rgba(76,201,160,0.15)',
            borderColor: 'transparent',
            pointRadius: 0,
            tension: 0.3,
          },
          {
            label: 'pc25',
            data: d.pc25,
            borderColor: 'rgba(76,201,160,0.35)',
            borderWidth: 1,
            pointRadius: 0,
            fill: false,
            tension: 0.3,
          },
          {
            label: 'Median',
            data: d.pc50,
            borderColor: '#4cc9a0',
            borderWidth: 2,
            pointRadius: 0,
            fill: false,
            tension: 0.3,
          },
          ...(d.kd_lo && d.kd_hi ? [{
            label: 'Kd −1σ',
            data: d.kd_lo,
            borderColor: 'rgba(200,180,80,0.45)',
            borderWidth: 1,
            borderDash: [3, 3],
            pointRadius: 0,
            fill: false,
            tension: 0.3,
          }, {
            label: 'Kd +1σ',
            data: d.kd_hi,
            borderColor: 'rgba(200,180,80,0.45)',
            borderWidth: 1,
            borderDash: [3, 3],
            pointRadius: 0,
            fill: false,
            tension: 0.3,
          }] : []),
        ]
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        interaction: { mode: 'index', intersect: false },
        plugins: {
          legend: { labels: { color: '#6e7f8d', font: { family: 'DM Mono', size: 11 } } },
          tooltip: { backgroundColor: '#161b22', borderColor: '#2a3441', borderWidth: 1,
                     titleColor: '#cdd9e5', bodyColor: '#6e7f8d',
                     titleFont: { family: 'DM Mono' }, bodyFont: { family: 'DM Mono' } }
        },
        scales: {
          x: { reverse: true, ticks: { color: '#6e7f8d', font: { family: 'DM Mono', size: 10 } },
               grid: { color: 'rgba(42,52,65,0.6)' },
               title: { display: true, text: 'Age (yrs BP)', color: '#6e7f8d', font: { family: 'DM Mono' } } },
          y: { ticks: { color: '#6e7f8d', font: { family: 'DM Mono', size: 10 } },
               grid: { color: 'rgba(42,52,65,0.6)' },
               title: { display: true,
                        text: d.mode==='semi' ? 'Relative drip rate (% of reference)' : 'Drip rate (min⁻¹)',
                        color: '#6e7f8d', font: { family: 'DM Mono' } } }
        }
      }
    });
  });
}

// ── Downloads ────────────────────────────────────────────────────────────────
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

function updateTheoKp(prefix) {
  const kpInput = document.getElementById(prefix + '_Kp');
  const theoField = document.getElementById(prefix + '_Kp_theo');
  const val = parseFloat(kpInput.value);
  if (val === -1 || kpInput.value.trim() === '-1') {
    // Get element for this TE
    const elem = document.getElementById(prefix + '_elem').value;
    if (THEO_KP[elem] !== undefined) {
      theoField.value = 'Kp = ' + THEO_KP[elem] + ' (Lindeman 2022)';
    } else {
      theoField.value = 'Kp = theoretical (Wang & Xu 2001)';
    }
    theoField.style.opacity = '0.7';
  } else {
    theoField.value = '';
    theoField.style.opacity = '0.4';
  }
}

// ── Labile fraction auto-calculation ────────────────────────────────────────────
function updateLabile(prefix) {
  const xF = parseFloat(document.getElementById(prefix + '_F').value) || 0;
  const xI = parseFloat(document.getElementById(prefix + '_InertF').value) || 0;
  const xL = Math.max(0, 1.0 - xI - xF);
  document.getElementById(prefix + '_labile').value = xL.toFixed(4);
  // Warn if fractions exceed 1
  const warn = (xI + xF) > 1.0;
  document.getElementById(prefix + '_F').style.borderColor    = warn ? 'var(--danger)' : '';
  document.getElementById(prefix + '_InertF').style.borderColor = warn ? 'var(--danger)' : '';
}

// ── Molecular weight lookup on element change ────────────────────────────────
const MOL_WT = {
  'Ni': 58.693, 'Co': 58.933, 'Cu': 63.546, 'V': 50.942,
  'Zn': 65.38,  'Cd': 112.411,'Al': 26.982, 'Pb': 207.2,
  'La': 138.905,'Ce': 140.116,'other': null
};
// Speleothem-specific inorganic Kd values from Lindeman et al. GCA 2022
const KP_DEFAULT = {
  'Ni': 1.1, 'Co': 4.4, 'Cu': 44,
  'V': -1, 'Zn': -1, 'Cd': -1, 'Al': -1,
  'Pb': -1, 'La': -1, 'Ce': -1, 'other': -1
};
function updateMolWt(prefix, elem) {
  const wt = MOL_WT[elem];
  if (wt !== null && wt !== undefined) {
    document.getElementById(prefix + '_mol_wt').value = wt;
  }
  // Set Kp default for this element
  const kp = KP_DEFAULT[elem] !== undefined ? KP_DEFAULT[elem] : -1;
  document.getElementById(prefix + '_Kp').value = kp;
  updateTheoKp(prefix);
}

// ── Utilities ────────────────────────────────────────────────────────────────
function fmtTime(s) {
  if (s < 60) return Math.round(s) + 's';
  return Math.round(s/60) + 'min';
}
function escHtml(s) {
  return s.replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
}
</script>
</body>
</html>
'''

@app.route('/')
def index():
    return Response(HTML, mimetype='text/html')


@app.route('/upload', methods=['POST'])
def upload():
    """Accept CSV uploads for depth/age, trace elements, and isotopes."""
    results = {}
    for key in ['depth_age', 'trace_elem1', 'trace_elem2', 'isotope1', 'isotope2']:
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
                # For depth_age and trace_elem1, return full data
                if key in ('depth_age', 'trace_elem1'):
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
            df = groups.mean(numeric_only=False).reset_index(drop=True)

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
        return 'File not found', 404
    return send_file(path, as_attachment=True)


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

            def _unit_to_ppb(unit):
                """Return multiplier to convert data in `unit` to ppb."""
                return {'ppb': 1.0, 'ppm': 1000.0,
                        'ug/g': 1000.0, 'mg/kg': 1000.0}.get(str(unit).strip(), 1.0)

            def load(key, col_depth, col_proxy, unit='ppb'):
                path = os.path.join(UPLOAD_FOLDER, key + '.csv')
                df = pd.read_csv(path)
                x = df[col_depth].to_numpy(dtype=float)
                y = df[col_proxy].to_numpy(dtype=float) * _unit_to_ppb(unit)
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
                    p25, p50, p75 = np.percentile(PD.proxy, [25, 50, 75])
                    limits = [0., p50 + 10. * (p75 - p25)]
                else:
                    limits = PDist.get_limits(DWF, PD, 3.)
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

        from drip_rate_util import driprates

        _analysis_mode = params.get('analysis_mode', 'full')
        log(f'Analysis mode: {_analysis_mode}')

        def _aq(row):
            """aq_conc converted to ppb — returns 1.0 in semi mode (shape only)."""
            if _analysis_mode == 'semi':
                return 1.0
            return float(params[row + '_aq_conc']) * _unit_to_ppb(params.get(row + '_aq_unit', 'ppb'))

        def make_te(pdist, row):
            return {
                'elem':    params[row + '_elem'],
                'mol_wt':  float(params[row + '_mol_wt']),
                'Kp':      float(params[row + '_Kp']),
                'Temp_C':  float(params['temp_C']),
                'F':       float(params[row + '_F']),
                'InertF':  float(params[row + '_InertF']),
                'aq_conc': np.float64(_aq(row)),
                'ca_conc': np.float64(
                    1.0 if _analysis_mode == 'semi'
                    else float(params['ca_conc']) * _unit_to_ppb(params.get('ca_unit', 'ppb'))
                ),
                'PDist':   copy.deepcopy(pdist),
            }

        # Loop over all TE PDists
        V_pdf_list = []
        for _i, _pd in enumerate(PDist_TEs):
            _rk = f'te{_i+1}'
            _pct = 65 + int(15 * _i / len(PDist_TEs))
            _set_stage(f'Computing drip rate PDF ({_rk.upper()})', _pct)
            log(f'Computing drip rate PDF for {_rk.upper()} …')
            _TE = make_te(_pd, _rk)
            _Ke = float(params.get(f'K_e{_i+1}', 1))
            _vpdf, V_age, V_span = driprates(
                float(params[f'{_rk}_Kd_mn']), float(params[f'{_rk}_Kd_sd']),
                _Ke, TE=_TE, calib=False)
            V_pdf_list.append(_vpdf)
        V_pdf_TE1 = V_pdf_list[0]
        V_pdf_TE2 = V_pdf_list[1] if len(V_pdf_list) > 1 else None

        # ── Kd sensitivity: run at Kd_mn ± Kd_sd (near-deterministic) ──
        _kd_sens_lo_pdfs = []
        _kd_sens_hi_pdfs = []
        for _i, _pd in enumerate(PDist_TEs):
            _rk = f'te{_i+1}'
            _TE_s = make_te(_pd, _rk)
            _Ke   = float(params.get(f'K_e{_i+1}', 1))
            _kmn  = float(params[f'{_rk}_Kd_mn'])
            _ksd  = float(params[f'{_rk}_Kd_sd'])
            _eps  = 0.001  # near-zero sd → near-deterministic
            _vlo, _, _ = driprates(_kmn - _ksd, _eps, _Ke, TE=_TE_s, calib=False)
            _vhi, _, _ = driprates(_kmn + _ksd, _eps, _Ke, TE=_TE_s, calib=False)
            _kd_sens_lo_pdfs.append(_vlo)
            _kd_sens_hi_pdfs.append(_vhi)
        log('Kd sensitivity passes complete')

        # ── Kd sensitivity joint PDFs ───────────────────────────────────────
        def _joint_and_median(pdf_list, rsw):
            n = len(pdf_list)
            if n == 1:
                jp = pdf_list[0].copy(); jp[jp < 0] = 0.
            else:
                jp = pdf_list[0].copy()
                for _v in pdf_list[1:]: jp = jp * _v
                jp[jp < 0] = 0.
                jp = np.power(jp, 1.0 / n)
            C = (jp.T * rsw).sum(axis=1)
            C[C == 0] = 1.
            jp = jp / C
            med = np.zeros(jp.shape[1])
            for _i in range(jp.shape[1]):
                cdf = np.cumsum(rsw * jp[:, _i]); cdf /= cdf[-1]
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
            V_pdf[V_pdf < 0.] = 0.
        else:
            log(f'Joint PDF: geometric mean of {_n_te_pdf} TEs')
            _prod = V_pdf_list[0].copy()
            for _v in V_pdf_list[1:]:
                _prod = _prod * _v
            _prod[_prod < 0.] = 0.
            V_pdf = np.power(_prod, 1.0 / _n_te_pdf)
        C = (V_pdf.T * V_rsw).sum(axis=1)
        V_pdf = V_pdf / C

        # ── 8. Summary statistics ───────────────────────────────────────────
        _set_stage('Extracting percentiles', 88)
        pcs = [5, 10, 25, 50, 75, 90, 95]
        V_pcs = {p: np.zeros(len(V_age)) for p in pcs}
        for i in range(V_pdf.shape[1]):
            cdf = np.cumsum(V_rsw * V_pdf[:, i])
            cdf /= cdf[-1]
            f_ = interp1d(cdf, V_span, kind='linear',
                          bounds_error=False, fill_value=(V_span[0], V_span[-1]))
            for p in pcs:
                V_pcs[p][i] = f_(p / 100.)

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

        # ── 10. Save outputs ────────────────────────────────────────────────
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
        chart_data = {
            'mode':  _analysis_mode,
            'age':   V_age.tolist(),
            'pc05':  V_pcs[5].tolist(),
            'pc25':  V_pcs[25].tolist(),
            'pc50':  V_pcs[50].tolist(),
            'pc75':  V_pcs[75].tolist(),
            'pc95':  V_pcs[95].tolist(),
            'kd_lo': _kd_med_lo.tolist(),
            'kd_hi': _kd_med_hi.tolist(),
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
