#!/usr/bin/env python3
"""
Quality Control for Sanger Chromatograms
Analyzes .ab1 files and generates QC report with chromatogram visualization
"""

import os
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.SeqIO import AbiIO
import pandas as pd
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import base64
from io import BytesIO
import webbrowser
import argparse

def print_header(text):
    """Print a formatted header"""
    width = 70
    print("\n" + "=" * width)
    print(f"  {text}")
    print("=" * width)

def print_step(step_num, total_steps, description):
    """Print a step indicator"""
    print(f"\n[Step {step_num}/{total_steps}] {description}")
    print("-" * 70)

def print_success(message):
    """Print a success message"""
    print(f"✓ {message}")

def print_info(message, indent=True):
    """Print an info message"""
    prefix = "  " if indent else ""
    print(f"{prefix}{message}")

def open_in_browser(file_path):
    """Open HTML file in default web browser (cross-platform)"""
    try:
        # Convert to absolute path and use file:// URL
        abs_path = Path(file_path).resolve()
        webbrowser.open(f'file://{abs_path}')
        return True
    except Exception as e:
        print(f"Could not auto-open browser: {e}")
        return False

def extract_trace_data(ab1_file):
    """Extract complete trace data for interactive viewer"""
    try:
        with open(ab1_file, 'rb') as f:
            record = AbiIO.AbiIterator(f).__next__()

        # Get base calls and quality
        seq_record = SeqIO.read(ab1_file, 'abi')
        sequence = str(seq_record.seq)
        quality = seq_record.letter_annotations.get('phred_quality', [])

        # Get peak positions
        peak_positions = record.annotations['abif_raw'].get('PLOC1', [])

        # Get trace data for all channels
        channels_data = {}
        channels = {'DATA9': 'G', 'DATA10': 'A', 'DATA11': 'T', 'DATA12': 'C'}

        for channel, base in channels.items():
            if channel in record.annotations['abif_raw']:
                # Downsample trace data for performance (every 2nd point)
                trace = record.annotations['abif_raw'][channel][::2]
                # Convert to list (handles both numpy arrays and tuples)
                channels_data[base] = list(trace)

        # Downsample peak positions too (divide by 2 since we downsampled traces)
        peak_positions_downsampled = [int(p/2) for p in peak_positions] if peak_positions else []

        return {
            'sequence': sequence,
            'quality': quality,
            'peak_positions': peak_positions_downsampled[:len(sequence)],  # Match sequence length
            'traces': channels_data,
            'trace_length': len(channels_data.get('A', []))
        }
    except Exception as e:
        print(f"  Warning: Could not extract trace data: {e}")
        return None

def plot_chromatogram(ab1_file, start_base=50, num_bases=150):
    """Generate chromatogram plot with sequence overlay

    Shows middle region by default (skipping poor quality start/end)
    """
    try:
        with open(ab1_file, 'rb') as f:
            record = AbiIO.AbiIterator(f).__next__()

        # Get trace data and base calls
        channels = {'DATA9': 'G', 'DATA10': 'A', 'DATA11': 'T', 'DATA12': 'C'}
        colors = {'G': 'black', 'A': 'green', 'T': 'red', 'C': 'blue'}

        # Get base calls and positions
        seq_record = SeqIO.read(ab1_file, 'abi')
        sequence = str(seq_record.seq)
        quality = seq_record.letter_annotations.get('phred_quality', [])

        # Get peak positions (where each base was called)
        peak_positions = record.annotations['abif_raw'].get('PLOC1', [])

        # Determine region to show
        end_base = min(start_base + num_bases, len(sequence))
        if start_base >= len(sequence):
            start_base = max(0, len(sequence) - num_bases)
            end_base = len(sequence)

        # Get trace region
        if peak_positions and len(peak_positions) > end_base:
            trace_start = peak_positions[start_base] if start_base < len(peak_positions) else 0
            trace_end = peak_positions[end_base-1] if end_base-1 < len(peak_positions) else len(record.annotations['abif_raw']['DATA9'])
        else:
            # Approximate if no peak positions
            trace_start = start_base * 10
            trace_end = end_base * 10

        fig, ax = plt.subplots(figsize=(18, 5))

        # Plot traces
        for channel, base in channels.items():
            if channel in record.annotations['abif_raw']:
                trace = record.annotations['abif_raw'][channel][trace_start:trace_end]
                ax.plot(trace, color=colors[base], label=base, linewidth=0.8, alpha=0.7)

        # Add base calls on top
        if peak_positions and len(peak_positions) > end_base:
            for i in range(start_base, end_base):
                if i < len(sequence) and i < len(peak_positions):
                    pos = peak_positions[i] - trace_start
                    base = sequence[i]

                    # Color code by quality
                    qual = quality[i] if i < len(quality) else 0
                    if qual >= 30:
                        qual_color = 'darkgreen'
                    elif qual >= 20:
                        qual_color = 'orange'
                    else:
                        qual_color = 'red'

                    ax.text(pos, ax.get_ylim()[1] * 0.95, base,
                           ha='center', va='top', fontsize=8,
                           color=qual_color, fontweight='bold')

        ax.set_xlabel('Trace Data Position (arbitrary units)', fontsize=10)
        ax.set_ylabel('Fluorescence Signal Intensity', fontsize=10)
        ax.set_title(f'{ab1_file.name} - Region showing bases {start_base}-{end_base}\n(Base call colors: green=Q≥30, orange=Q≥20, red=Q<20)',
                    fontsize=11, fontweight='bold')
        ax.legend(loc='upper right')
        ax.grid(True, alpha=0.2)

        # Convert to base64
        buf = BytesIO()
        plt.savefig(buf, format='png', dpi=120, bbox_inches='tight')
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode('utf-8')
        plt.close()

        return img_base64
    except Exception as e:
        print(f"  Warning: Could not generate chromatogram plot: {e}")
        return None

def analyze_chromatogram(ab1_file, output_dir=None):
    """Analyze a single .ab1 chromatogram file"""
    try:
        record = SeqIO.read(ab1_file, "abi")

        # Extract quality scores
        quality_scores = record.letter_annotations.get("phred_quality", [])

        # Calculate QC metrics
        seq_length = len(record.seq)
        avg_quality = sum(quality_scores) / len(quality_scores) if quality_scores else 0
        high_quality_bases = sum(1 for q in quality_scores if q >= 20)
        low_quality_bases = sum(1 for q in quality_scores if q < 20)

        # Determine pass/fail
        qc_pass = avg_quality >= 20 and high_quality_bases / seq_length >= 0.8

        result = {
            "file": ab1_file.name,
            "length": seq_length,
            "avg_quality": round(avg_quality, 2),
            "high_quality_bases": high_quality_bases,
            "low_quality_bases": low_quality_bases,
            "percent_high_quality": round(100 * high_quality_bases / seq_length, 2),
            "qc_status": "PASS" if qc_pass else "FAIL",
            "sequence": str(record.seq),
            "quality_scores": quality_scores
        }

        # Generate chromatogram plot if output_dir provided
        if output_dir:
            result["chromatogram_img"] = plot_chromatogram(ab1_file)
            # Also extract full trace data for interactive viewer
            result["trace_data"] = extract_trace_data(ab1_file)

        return result
    except Exception as e:
        return {
            "file": ab1_file.name,
            "error": str(e),
            "qc_status": "ERROR"
        }

def generate_quality_heatmap(quality_scores):
    """Generate HTML for quality heatmap visualization"""
    if not quality_scores:
        return ""

    # Divide sequence into segments for heatmap
    segment_size = max(1, len(quality_scores) // 50)  # Show ~50 segments
    segments = []

    for i in range(0, len(quality_scores), segment_size):
        segment_quals = quality_scores[i:i+segment_size]
        avg_qual = sum(segment_quals) / len(segment_quals) if segment_quals else 0

        if avg_qual >= 30:
            qual_class = "q-high"
        elif avg_qual >= 20:
            qual_class = "q-medium"
        elif avg_qual > 0:
            qual_class = "q-low"
        else:
            qual_class = "q-none"

        segments.append(qual_class)

    heatmap_html = '<div class="quality-heatmap">'
    for qual_class in segments:
        heatmap_html += f'<div class="quality-segment {qual_class}"></div>'
    heatmap_html += '</div>'

    return heatmap_html

def generate_html_report(results, output_file):
    """Generate HTML report with chromatogram visualizations using new design system"""
    from datetime import datetime
    from pathlib import Path

    # Calculate summary statistics
    total = len(results)
    passed = sum(1 for r in results if r['qc_status'] == "PASS")
    failed = sum(1 for r in results if r['qc_status'] == "FAIL")
    errors = sum(1 for r in results if r['qc_status'] == "ERROR")

    # Average quality of passed sequences
    passed_results = [r for r in results if r['qc_status'] == "PASS" and 'avg_quality' in r]
    avg_quality = sum(r['avg_quality'] for r in passed_results) / len(passed_results) if passed_results else 0

    # Read CSS files and embed them
    project_root = Path(__file__).parent.parent.parent
    base_css = (project_root / "tracking/styles/base.css").read_text()
    components_css = (project_root / "tracking/styles/components.css").read_text()
    reports_css = (project_root / "tracking/styles/reports.css").read_text()
    chromatogram_css = (project_root / "tracking/styles/chromatogram.css").read_text()

    # Build HTML
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Quality Control Report - DNA Barcoding</title>

    <!-- Embedded CSS for reliable loading -->
    <style>
{base_css}

{components_css}

{reports_css}

{chromatogram_css}
    </style>

    <!-- Interactive Chromatogram Viewer JavaScript -->
    <script>
        // Render chromatogram on canvas
        function renderChromatogram(sampleId, startBase) {{
            const traceData = window['traceData_' + sampleId];
            if (!traceData) return;

            const canvas = document.getElementById('chromato-canvas-' + sampleId);
            if (!canvas) return;

            const ctx = canvas.getContext('2d');
            const width = canvas.width;
            const height = canvas.height;

            // Clear canvas
            ctx.clearRect(0, 0, width, height);

            // Configuration - Reduced window for better peak visibility
            const basesToShow = 60;  // Changed from 150 to 60 for clearer peaks
            const endBase = Math.min(startBase + basesToShow, traceData.sequence.length);
            const actualStart = Math.max(0, startBase);

            // Get trace window
            const peakPos = traceData.peak_positions;
            if (!peakPos || peakPos.length === 0) return;

            const traceStart = peakPos[actualStart] || 0;
            const traceEnd = peakPos[Math.min(endBase - 1, peakPos.length - 1)] || traceData.trace_length;

            // Draw traces
            const traces = {{
                'G': {{ color: 'rgba(0, 0, 0, 0.7)', data: traceData.traces.G }},
                'A': {{ color: 'rgba(0, 200, 0, 0.7)', data: traceData.traces.A }},
                'T': {{ color: 'rgba(255, 0, 0, 0.7)', data: traceData.traces.T }},
                'C': {{ color: 'rgba(0, 0, 255, 0.7)', data: traceData.traces.C }}
            }};

            // Find max intensity in window for scaling
            let maxIntensity = 0;
            for (const base in traces) {{
                const slice = traces[base].data.slice(traceStart, traceEnd);
                maxIntensity = Math.max(maxIntensity, ...slice);
            }}

            // Scale factor
            const yScale = (height * 0.75) / maxIntensity;
            const xScale = width / (traceEnd - traceStart);

            // Draw each trace with SMOOTH curves
            for (const base in traces) {{
                const trace = traces[base];
                ctx.strokeStyle = trace.color;
                ctx.lineWidth = 1.5;
                ctx.lineJoin = 'round';  // Smooth corners
                ctx.lineCap = 'round';   // Smooth line ends
                ctx.beginPath();

                // Convert trace data to points
                const points = [];
                for (let i = traceStart; i < traceEnd; i++) {{
                    const x = (i - traceStart) * xScale;
                    const y = height - (trace.data[i] * yScale);
                    points.push({{ x, y }});
                }}

                if (points.length > 0) {{
                    ctx.moveTo(points[0].x, points[0].y);

                    // Draw smooth curves using quadratic interpolation
                    for (let i = 1; i < points.length - 1; i++) {{
                        const xc = (points[i].x + points[i + 1].x) / 2;
                        const yc = (points[i].y + points[i + 1].y) / 2;
                        ctx.quadraticCurveTo(points[i].x, points[i].y, xc, yc);
                    }}

                    // Draw last segment
                    if (points.length > 1) {{
                        const last = points[points.length - 1];
                        const prev = points[points.length - 2];
                        ctx.quadraticCurveTo(prev.x, prev.y, last.x, last.y);
                    }}
                }}

                ctx.stroke();
            }}

            // Draw base calls
            ctx.font = 'bold 14px monospace';
            ctx.textAlign = 'center';

            for (let i = actualStart; i < endBase; i++) {{
                if (i >= traceData.sequence.length) break;

                const base = traceData.sequence[i];
                const qual = traceData.quality[i] || 0;

                // Position on canvas
                const peakIdx = peakPos[i];
                if (peakIdx === undefined || peakIdx < traceStart || peakIdx >= traceEnd) continue;

                const x = (peakIdx - traceStart) * xScale;
                const y = 20;

                // Color by quality
                if (qual >= 30) {{
                    ctx.fillStyle = '#006400';  // Dark green
                }} else if (qual >= 20) {{
                    ctx.fillStyle = '#ff8c00';  // Orange
                }} else {{
                    ctx.fillStyle = '#dc3545';  // Red
                }}

                ctx.fillText(base, x, y);
            }}

            // Update position display (minimal format)
            const posDisplay = document.getElementById('chromato-position-' + sampleId);
            if (posDisplay) {{
                const totalBases = traceData.sequence.length;
                posDisplay.innerHTML = `
                    <span style="color: var(--purple);">📍 Base:</span>
                    <span style="color: var(--text-primary); font-weight: 700;">${{actualStart}}-${{endBase}}</span>
                    <span style="color: var(--text-secondary);"> (of ${{totalBases}})</span>
                    <span style="margin: 0 0.5rem; color: var(--text-secondary);">•</span>
                    <span style="color: var(--cyan); font-weight: 600;">Trace:</span>
                    <span style="color: var(--text-primary); font-weight: 700;">${{traceStart}}-${{traceEnd}}</span>
                `;
            }}
        }}

        // Update chromatogram from slider
        function updateChromatogram(sampleId, value) {{
            const startBase = parseInt(value);
            renderChromatogram(sampleId, startBase);
        }}

        // Move chromatogram by offset
        function moveChromatogram(sampleId, offset) {{
            const slider = document.getElementById('chromato-slider-' + sampleId);
            const newValue = parseInt(slider.value) + offset;
            const maxValue = parseInt(slider.max);

            slider.value = Math.max(0, Math.min(newValue, maxValue));
            updateChromatogram(sampleId, slider.value);
        }}

        // Reset to default view
        function resetChromatogram(sampleId) {{
            const slider = document.getElementById('chromato-slider-' + sampleId);
            slider.value = 50;
            updateChromatogram(sampleId, 50);
        }}

        // Add CSS for controls
        const style = document.createElement('style');
        style.textContent = `
            .chromato-controls {{
                margin-top: 1rem;
                display: flex;
                flex-direction: column;
                gap: 1rem;
                padding: 1rem;
                background: var(--bg-secondary);
                border-radius: var(--radius-md);
            }}

            .control-group {{
                display: flex;
                gap: 0.5rem;
                flex-wrap: wrap;
            }}

            .slider-group {{
                display: flex;
                align-items: center;
                gap: 1rem;
                flex-wrap: wrap;
            }}

            .chromato-slider {{
                -webkit-appearance: none;
                appearance: none;
                height: 12px;
                background: linear-gradient(90deg, var(--purple), var(--pink));
                border-radius: 6px;
                outline: none;
                transition: opacity 0.2s;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            }}

            .chromato-slider:hover {{
                opacity: 0.9;
                box-shadow: 0 2px 8px rgba(189, 147, 249, 0.4);
            }}

            /* Chrome/Safari/Edge slider thumb - BIG and OBVIOUS! */
            .chromato-slider::-webkit-slider-thumb {{
                -webkit-appearance: none;
                appearance: none;
                width: 32px;
                height: 32px;
                background: linear-gradient(135deg, var(--purple) 0%, var(--pink) 100%);
                border: 4px solid white;
                border-radius: 50%;
                cursor: grab;
                box-shadow: 0 3px 8px rgba(0,0,0,0.3), 0 0 0 4px rgba(189, 147, 249, 0.2);
                transition: all 0.2s ease;
                animation: pulse 2s ease-in-out infinite;
            }}

            .chromato-slider::-webkit-slider-thumb:hover {{
                transform: scale(1.15);
                box-shadow: 0 4px 12px rgba(189, 147, 249, 0.6), 0 0 0 6px rgba(189, 147, 249, 0.3);
            }}

            .chromato-slider::-webkit-slider-thumb:active {{
                cursor: grabbing;
                background: linear-gradient(135deg, var(--pink) 0%, var(--purple) 100%);
                transform: scale(1.05);
                animation: none;
            }}

            /* Firefox slider thumb - BIG and OBVIOUS! */
            .chromato-slider::-moz-range-thumb {{
                width: 32px;
                height: 32px;
                background: linear-gradient(135deg, var(--purple) 0%, var(--pink) 100%);
                border: 4px solid white;
                border-radius: 50%;
                cursor: grab;
                box-shadow: 0 3px 8px rgba(0,0,0,0.3), 0 0 0 4px rgba(189, 147, 249, 0.2);
                transition: all 0.2s ease;
            }}

            .chromato-slider::-moz-range-thumb:hover {{
                transform: scale(1.15);
                box-shadow: 0 4px 12px rgba(189, 147, 249, 0.6), 0 0 0 6px rgba(189, 147, 249, 0.3);
            }}

            .chromato-slider::-moz-range-thumb:active {{
                cursor: grabbing;
                background: linear-gradient(135deg, var(--pink) 0%, var(--purple) 100%);
                transform: scale(1.05);
            }}

            /* Pulsing animation to draw attention */
            @keyframes pulse {{
                0%, 100% {{
                    box-shadow: 0 3px 8px rgba(0,0,0,0.3), 0 0 0 4px rgba(189, 147, 249, 0.2);
                }}
                50% {{
                    box-shadow: 0 3px 8px rgba(0,0,0,0.3), 0 0 0 8px rgba(189, 147, 249, 0.4);
                }}
            }}

            @media (max-width: 768px) {{
                .chromato-controls {{
                    padding: 0.75rem;
                }}

                .control-group {{
                    flex-direction: column;
                }}

                .slider-group {{
                    flex-direction: column;
                    align-items: flex-start;
                }}

                .chromato-slider {{
                    width: 100% !important;
                    max-width: none !important;
                }}
            }}
        `;
        document.head.appendChild(style);
    </script>
</head>
<body>
    <!-- Report Header -->
    <header class="report-header">
        <h1>🔬 Quality Control Report</h1>
        <div class="progress-badge">Step 1 of 5</div>
        <div class="report-date">Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</div>
    </header>

    <!-- Summary Dashboard -->
    <section class="summary-dashboard">
        <div class="metric-card metric-success">
            <div class="metric-value">{passed}</div>
            <div class="metric-label">Passed QC</div>
            <div class="metric-icon">✓</div>
        </div>

        <div class="metric-card metric-{"fail" if failed > 0 else "success"}">
            <div class="metric-value">{failed}</div>
            <div class="metric-label">Failed QC</div>
            <div class="metric-icon">{"✗" if failed > 0 else "✓"}</div>
        </div>

        <div class="metric-card metric-primary">
            <div class="metric-value">{avg_quality:.1f}</div>
            <div class="metric-label">Avg. Quality (Q)</div>
            <div class="metric-icon">📊</div>
        </div>

        <div class="metric-card metric-{"fail" if errors > 0 else "success"}">
            <div class="metric-value">{errors}</div>
            <div class="metric-label">Errors</div>
            <div class="metric-icon">{"⚠" if errors > 0 else "✓"}</div>
        </div>
    </section>

    <!-- Main Content -->
    <main class="report-content">
        <div class="content-section">
            <h2>Quality Control Results</h2>
            <p>Sanger sequencing chromatogram analysis with Phred quality scores</p>

            <div class="info-box info-tip">
                <strong>💡 Quality Thresholds:</strong>
                <ul>
                    <li><strong>Q ≥ 30:</strong> High quality (99.9% base call accuracy)</li>
                    <li><strong>Q ≥ 20:</strong> Acceptable quality (99% accuracy)</li>
                    <li><strong>Q &lt; 20:</strong> Low quality (may need review)</li>
                    <li><strong>Pass criteria:</strong> ≥80% bases with Q≥20, average Q≥20</li>
                </ul>
            </div>

            <div class="info-box info-success" style="margin-top: 1rem; background: linear-gradient(135deg, rgba(189, 147, 249, 0.1), rgba(255, 121, 198, 0.1));">
                <strong>🔬 Interactive Chromatogram Viewer Guide:</strong>
                <ul style="margin: 0.5rem 0 0 1.5rem;">
                    <li><strong>👆 Drag the BIG circular slider</strong> to scroll through the entire sequence</li>
                    <li><strong>⏮️ ⏭️ Use arrow buttons</strong> to jump forward/backward by 20 bases</li>
                    <li><strong>🎨 Trace colors:</strong> A=green, T=red, C=blue, G=black</li>
                    <li><strong>📊 Base quality on peaks:</strong> Green text=high (Q≥30), Orange=medium (Q≥20), Red=low (Q&lt;20)</li>
                </ul>
                <div style="margin-top: 1rem; padding: 0.75rem; background: var(--bg-cyan-light); border-radius: var(--radius-md); border-left: 3px solid var(--cyan);">
                    <strong>📍 Understanding Positions:</strong>
                    <ul style="margin: 0.5rem 0 0 1.5rem;">
                        <li><strong>Base Position:</strong> Which DNA bases (A, T, G, C) you're viewing (e.g., bases 50-110)</li>
                        <li><strong>Trace Position:</strong> The raw fluorescence data points from the sequencer (each base has ~10-15 trace points)</li>
                    </ul>
                </div>
            </div>
"""

    # Generate results for each sequence
    for result in results:
        filename = result['file']
        status = result['qc_status']

        if status == "ERROR":
            html += f"""
            <div class="collapsible-section">
                <button class="collapsible-toggle">
                    <span class="toggle-icon">▶</span>
                    <strong>{filename}</strong>
                    <span class="badge badge-fail">ERROR</span>
                </button>
                <div class="collapsible-content">
                    <div class="info-box info-error">
                        <strong>⚠️ Processing Error:</strong> {result.get('error', 'Unknown error')}
                    </div>
                </div>
            </div>
"""
            continue

        # Get metrics
        length = result.get('length', 0)
        avg_qual = result.get('avg_quality', 0)
        high_qual_pct = result.get('percent_high_quality', 0)

        # Determine badge and styling
        if status == "PASS":
            badge_class = "badge-success"
            badge_text = "PASS"
            status_class = "pass"
            status_upper = "PASS"
        else:
            badge_class = "badge-fail"
            badge_text = "FAIL"
            status_class = "fail"
            status_upper = "FAIL"

        html += f"""
            <div class="collapsible-section">
                <button class="collapsible-toggle">
                    <span class="toggle-icon">▶</span>
                    <strong>{filename}</strong>
                    <span class="badge {badge_class}">{badge_text}</span>
                    <em style="margin-left: auto; color: var(--text-secondary);">{length} bp | Q{avg_qual:.1f}</em>
                </button>
                <div class="collapsible-content">
                    <!-- Quality Metrics Summary Table -->
                    <div style="overflow-x: auto; margin-bottom: 1.5rem;">
                        <table class="data-table" style="margin: 0;">
                            <thead>
                                <tr style="background: linear-gradient(135deg, var(--bg-green-light), var(--bg-cyan-light));">
                                    <th>📏 Length</th>
                                    <th>⭐ Avg Quality</th>
                                    <th>✓ High Quality %</th>
                                    <th>📊 Q≥20 Bases</th>
                                    <th>📉 Q<20 Bases</th>
                                    <th>🎯 Status</th>
                                </tr>
                            </thead>
                            <tbody>
                                <tr class="row-{status_class}">
                                    <td><strong>{length} bp</strong></td>
                                    <td><span class="badge badge-{"success" if avg_qual >= 30 else "warning" if avg_qual >= 20 else "fail"}">Q {avg_qual:.1f}</span></td>
                                    <td><strong>{high_qual_pct:.1f}%</strong></td>
                                    <td>{result.get('high_quality_bases', 0)} / {length}</td>
                                    <td>{result.get('low_quality_bases', 0)} / {length}</td>
                                    <td><span class="badge badge-{status_class}">{status_upper}</span></td>
                                </tr>
                            </tbody>
                        </table>
                    </div>

                    <!-- Quality Overview Heatmap -->
                    <div class="quality-overview">
                        <div class="quality-overview-label">Quality Distribution Across Sequence</div>
                        {generate_quality_heatmap(result.get('quality_scores', []))}
                    </div>
"""

        # Add interactive chromatogram if trace data available
        if "trace_data" in result and result["trace_data"]:
            trace_data = result["trace_data"]
            trace_data_json = json.dumps(trace_data)

            html += f"""
                    <!-- Interactive Chromatogram Visualization -->
                    <div class="chromatogram-wrapper">
                        <div class="chromatogram-header">
                            <h4 class="chromatogram-title">🔬 Interactive Chromatogram Viewer</h4>
                            <div class="chromatogram-meta">
                                <span class="badge badge-success">NEW!</span>
                                <span style="margin-left: 1rem;">Drag slider to explore • Full sequence: {len(trace_data['sequence'])} bases</span>
                            </div>
                        </div>

                        <div class="chromatogram-viewer">
                            <!-- Canvas for rendering traces -->
                            <canvas id="chromato-canvas-{filename.replace('.', '_')}"
                                    width="1400"
                                    height="350"
                                    style="width: 100%; border: 1px solid var(--border-color); border-radius: var(--radius-md); background: white; cursor: crosshair;">
                            </canvas>

                            <!-- Navigation Controls -->
                            <div class="chromato-controls">
                                <div class="control-group">
                                    <button class="btn btn-nav" onclick="moveChromatogram('{filename.replace('.', '_')}', -20)">
                                        ⏮️ Back 20
                                    </button>
                                    <button class="btn btn-nav" onclick="moveChromatogram('{filename.replace('.', '_')}', 20)">
                                        Forward 20 ⏭️
                                    </button>
                                    <button class="btn btn-secondary" onclick="resetChromatogram('{filename.replace('.', '_')}')">
                                        🔄 Reset
                                    </button>
                                </div>

                                <div class="slider-group">
                                    <label for="chromato-slider-{filename.replace('.', '_')}" style="font-weight: 600; margin-right: 1rem; font-size: 0.95rem; display: flex; align-items: center; gap: 0.5rem; flex-wrap: wrap;">
                                        <span id="chromato-position-{filename.replace('.', '_')}" style="font-size: 0.95rem;">
                                            <span style="color: var(--purple);">📍 Base:</span>
                                            <span style="color: var(--text-primary); font-weight: 700;">50-110</span>
                                            <span style="color: var(--text-secondary);"> (of {len(trace_data['sequence'])})</span>
                                            <span style="margin: 0 0.5rem; color: var(--text-secondary);">•</span>
                                            <span style="color: var(--cyan); font-weight: 600;">Trace:</span>
                                            <span style="color: var(--text-primary); font-weight: 700;">Loading...</span>
                                        </span>
                                    </label>
                                    <input type="range"
                                           id="chromato-slider-{filename.replace('.', '_')}"
                                           class="chromato-slider"
                                           min="0"
                                           max="{max(0, len(trace_data['sequence']) - 60)}"
                                           value="50"
                                           oninput="updateChromatogram('{filename.replace('.', '_')}', this.value)"
                                           style="width: 100%; max-width: 600px;">
                                </div>
                            </div>
                        </div>

                        <!-- Trace Legend - Inline -->
                        <div style="margin-top: 1rem; padding: 0.75rem; background: var(--bg-secondary); border-radius: var(--radius-md); text-align: center; font-size: 0.9rem;">
                            <strong style="color: var(--text-secondary);">Traces:</strong>
                            <span style="margin: 0 1rem; color: var(--green); font-weight: 600;">A</span>
                            <span style="margin: 0 1rem; color: var(--red); font-weight: 600;">T</span>
                            <span style="margin: 0 1rem; color: var(--cyan); font-weight: 600;">C</span>
                            <span style="margin: 0 1rem; color: black; font-weight: 600;">G</span>
                        </div>
                    </div>

                    <!-- JavaScript for this chromatogram -->
                    <script>
                        // Store trace data for this sample
                        window.traceData_{filename.replace('.', '_')} = {trace_data_json};

                        // Initialize the chromatogram viewer
                        document.addEventListener('DOMContentLoaded', function() {{
                            renderChromatogram('{filename.replace('.', '_')}', 50);
                        }});
                    </script>
"""
        # Fallback to static image if trace data not available
        elif "chromatogram_img" in result and result["chromatogram_img"]:
            html += f"""
                    <!-- Static Chromatogram (trace data unavailable) -->
                    <div class="chromatogram-wrapper">
                        <div class="chromatogram-header">
                            <h4 class="chromatogram-title">Chromatogram View</h4>
                            <div class="chromatogram-meta">Showing bases 50-200 (middle region)</div>
                        </div>
                        <div class="chromatogram-viewer">
                            <div class="chromatogram-canvas-container">
                                <img src="data:image/png;base64,{result['chromatogram_img']}"
                                     alt="Chromatogram for {filename}"
                                     class="chromatogram-image">
                            </div>
                        </div>
                    </div>
"""

        # Add full sequence
        sequence = result.get('sequence', 'N/A')
        if sequence and sequence != 'N/A':
            html += f"""
                    <!-- Full Sequence -->
                    <h4 style="margin-top: 1.5rem;">Full Sequence ({len(sequence)} bp)</h4>
                    <div class="sequence-display">{sequence}</div>
"""

        html += """
                </div>
            </div>
"""

    html += """
        </div>
    </main>

    <!-- Footer with Help -->
    <footer class="report-footer">
        <div class="help-section">
            <h3>Understanding Quality Control</h3>

            <h4>What are Phred Quality Scores?</h4>
            <p>Phred quality scores indicate the confidence in each base call. Q20 means 99% accuracy (1 in 100 chance of error), Q30 means 99.9% accuracy (1 in 1000 chance of error).</p>

            <h4>Why do some sequences fail?</h4>
            <ul>
                <li><strong>Too short:</strong> Sequences under 500bp may not provide enough information for species identification</li>
                <li><strong>Low quality:</strong> Poor chromatogram peaks can lead to incorrect base calls</li>
                <li><strong>Contamination:</strong> Mixed signals from multiple DNA sources</li>
                <li><strong>Degraded DNA:</strong> Old or improperly stored samples</li>
            </ul>

            <h4>How to improve sequence quality</h4>
            <ol>
                <li>Ensure DNA is fresh and properly stored (-20°C or -80°C)</li>
                <li>Use high-quality primers designed for your target region</li>
                <li>Optimize PCR conditions (annealing temperature, cycle number)</li>
                <li>Clean PCR products before sequencing (remove primers and dNTPs)</li>
                <li>Request re-sequencing for failed samples</li>
            </ol>

            <h4>What does the chromatogram show?</h4>
            <p>The chromatogram displays raw fluorescence signals from the sequencer. Clean, sharp peaks with good separation indicate high-quality sequence data. Overlapping peaks or noisy baselines suggest quality issues.</p>
        </div>
    </footer>

    <!-- JavaScript for collapsible sections -->
    <script>
        document.querySelectorAll('.collapsible-toggle').forEach(button => {
            button.addEventListener('click', () => {
                const section = button.parentElement;
                section.classList.toggle('expanded');
            });
        });
    </script>
</body>
</html>
"""

    with open(output_file, 'w') as f:
        f.write(html)

def main():
    """Main QC function"""
    # Parse arguments
    parser = argparse.ArgumentParser(
        description='Quality Control for Sanger Chromatograms',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python qc_chromatograms.py data/my_sequences/
  python qc_chromatograms.py data/my_sequences/ results/ --open
        """
    )
    parser.add_argument('input_dir', type=str, help='Directory containing .ab1 files')
    parser.add_argument('output_dir', type=str, nargs='?', default='results',
                       help='Output directory for results (default: results/)')
    parser.add_argument('--open', action='store_true',
                       help='Automatically open HTML report in web browser')

    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Header
    print_header("DNA BARCODING QUALITY CONTROL")
    print_info(f"Analyzing Sanger chromatograms (.ab1 files)", indent=False)
    print_info(f"Input: {input_dir}", indent=False)
    print_info(f"Output: {output_dir}", indent=False)

    # Step 1: Find files
    print_step(1, 3, "Finding chromatogram files")
    ab1_files = list(input_dir.glob("*.ab1"))

    if not ab1_files:
        print(f"\nERROR: No .ab1 files found in {input_dir}")
        print("Please check that your chromatogram files are in the correct directory.")
        sys.exit(1)

    # Sort files to group forward/reverse pairs together
    def sort_key(filepath):
        stem = filepath.stem
        parts = stem.rsplit('_', 1)
        if len(parts) == 2:
            sample_name, direction = parts
            return (sample_name, direction)
        return (stem, '')

    ab1_files = sorted(ab1_files, key=sort_key)
    print_success(f"Found {len(ab1_files)} chromatogram files")
    for f in ab1_files:
        print_info(f"• {f.name}")

    # Step 2: Analyze each file
    print_step(2, 3, f"Analyzing chromatograms")
    print_info("Checking sequence quality, length, and reading frames...")

    results = []
    passed_sequences = []

    for i, ab1_file in enumerate(ab1_files, 1):
        print_info(f"[{i}/{len(ab1_files)}] {ab1_file.name}...", indent=True)
        result = analyze_chromatogram(ab1_file, output_dir=output_dir)
        results.append(result)

        # Save passed sequences to FASTA
        if result['qc_status'] == "PASS" and 'sequence' in result:
            passed_sequences.append({
                'id': ab1_file.stem,
                'seq': result['sequence']
            })
            print_info(f"    ✓ PASSED (length: {result['length']} bp, avg quality: {result['avg_quality']:.1f})", indent=True)
        else:
            print_info(f"    ✗ FAILED ({result.get('fail_reason', 'Unknown reason')})", indent=True)

    # Step 3: Generate outputs
    print_step(3, 3, "Generating reports")

    # HTML report
    html_file = output_dir / "qc_report.html"
    print_info("Creating HTML report with chromatogram visualizations...")
    generate_html_report(results, html_file)
    print_success(f"HTML report: {html_file}")

    # JSON report
    json_file = output_dir / "qc_results.json"
    with open(json_file, 'w') as f:
        json.dump(results, f, indent=2)
    print_success(f"JSON data: {json_file}")

    # FASTA file for passed sequences
    if passed_sequences:
        fasta_file = output_dir / "passed_sequences.fasta"
        with open(fasta_file, 'w') as f:
            for seq in passed_sequences:
                f.write(f">{seq['id']}\n{seq['seq']}\n")
        print_success(f"Passed sequences: {fasta_file}")

    # Final Summary
    passed = sum(1 for r in results if r['qc_status'] == "PASS")
    failed = sum(1 for r in results if r['qc_status'] == "FAIL")
    errors = sum(1 for r in results if r['qc_status'] == "ERROR")

    print_header("QUALITY CONTROL COMPLETE")
    print_info(f"Total sequences analyzed: {len(results)}", indent=False)
    print_info(f"✓ Passed QC: {passed} sequences", indent=False)
    print_info(f"✗ Failed QC: {failed} sequences", indent=False)
    if errors > 0:
        print_info(f"⚠ Errors: {errors} sequences", indent=False)

    print("\n" + "=" * 70)
    print("  NEXT STEPS:")
    print("=" * 70)
    print_info("1. Open the HTML report to view detailed results:", indent=False)
    print_info(f"   {html_file.resolve()}", indent=False)
    print_info("", indent=False)
    print_info("2. Passed sequences are ready for alignment:", indent=False)
    if passed_sequences:
        print_info(f"   {fasta_file.resolve()}", indent=False)
    else:
        print_info("   (No sequences passed QC)", indent=False)
    print("=" * 70 + "\n")

    # Auto-open browser if requested
    if args.open:
        print_info("Opening HTML report in your web browser...", indent=False)
        if open_in_browser(html_file):
            print_success("Report opened successfully")
        else:
            print_info(f"Please open manually: {html_file.resolve()}", indent=False)

if __name__ == "__main__":
    main()
