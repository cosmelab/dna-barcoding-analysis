#!/usr/bin/env python3
"""
Run Lab Data Analysis - Wrapper script.

Runs all analysis scripts in order and generates the combined HTML report.
Output: results/lab_analysis/lab_analysis_report.html
"""

import sys
from pathlib import Path

# Add module directory to path
module_dir = Path(__file__).parent
sys.path.insert(0, str(module_dir))


def main():
    """Run all analysis scripts."""
    print()
    print("=" * 70)
    print("LAB DATA ANALYSIS - ENTM 201L DNA BARCODING")
    print("=" * 70)
    print()

    # 1. Extraction comparison
    print("[1/7] Extraction comparison...")
    from plot_extraction import main as extraction_main
    extraction_main()
    print()

    # 2. Quality assessment
    print("[2/7] DNA quality assessment...")
    from plot_quality import main as quality_main
    quality_main()
    print()

    # 3. PCR results
    print("[3/7] PCR results...")
    from plot_pcr import main as pcr_main
    pcr_main()
    print()

    # 4. Sequencing results
    print("[4/7] Sequencing results...")
    from plot_sequencing import main as sequencing_main
    sequencing_main()
    print()

    # 5. Team comparison
    print("[5/7] Team comparison...")
    from plot_teams import main as teams_main
    teams_main()
    print()

    # 6. Pipeline success correlation
    print("[6/7] Pipeline success correlation...")
    from plot_pipeline_success import main as pipeline_main
    pipeline_main()
    print()

    # 7. Generate combined report
    print("[7/7] Generating HTML report...")
    from generate_report import main as report_main
    report_main()
    print()

    # Summary
    output_dir = Path(__file__).parent.parent.parent / "results" / "lab_analysis"

    print()
    print("=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print()
    print(f"Output directory: {output_dir}")
    print()
    print("Generated files:")
    for f in sorted(output_dir.glob("*")):
        if f.is_file():
            size = f.stat().st_size
            if size > 1024 * 1024:
                size_str = f"{size / (1024*1024):.1f} MB"
            elif size > 1024:
                size_str = f"{size / 1024:.1f} KB"
            else:
                size_str = f"{size} B"
            print(f"  {f.name} ({size_str})")
    print()
    print("Open the report:")
    print(f"  {output_dir / 'lab_analysis_report.html'}")
    print()


if __name__ == "__main__":
    main()
