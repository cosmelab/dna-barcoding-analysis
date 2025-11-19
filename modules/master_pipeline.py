#!/usr/bin/env python3
"""
DNA Barcoding Master Pipeline
Runs complete analysis: QC -> Alignment -> Phylogeny -> Identification
"""

import os
import sys
import subprocess
from pathlib import Path
from datetime import datetime

def run_command(cmd, step_name):
    """Run a command and handle errors"""
    print(f"\n{'='*60}")
    print(f"Step: {step_name}")
    print(f"{'='*60}")
    print(f"Running: {' '.join(cmd)}\n")

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(result.stdout)
        if result.stderr:
            print("Warnings:", result.stderr)
        print(f"✓ {step_name} completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"✗ {step_name} failed")
        print("Error output:", e.stderr)
        return False

def main():
    """Run complete DNA barcoding pipeline"""
    if len(sys.argv) < 2:
        print("Usage: python master_pipeline.py <input_directory> [output_directory]")
        print("Example: python master_pipeline.py data/my_sequences/ results/")
        print("\nThis will run:")
        print("  1. Quality Control")
        print("  2. Species Identification (BLAST)")
        sys.exit(1)

    input_dir = Path(sys.argv[1])
    output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("results")

    # Create output directory
    output_dir.mkdir(exist_ok=True)

    # Create timestamped run directory
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = output_dir / f"run_{timestamp}"
    run_dir.mkdir(exist_ok=True)

    print(f"\nDNA Barcoding Analysis Pipeline")
    print(f"{'='*60}")
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {run_dir}")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Step 1: Quality Control
    qc_script = Path(__file__).parent / "01_quality_control" / "qc_chromatograms.py"
    qc_success = run_command(
        ["python", str(qc_script), str(input_dir), str(run_dir)],
        "Quality Control"
    )

    if not qc_success:
        print("\n✗ Pipeline failed at Quality Control step")
        sys.exit(1)

    # Check if we have passed sequences
    passed_fasta = run_dir / "passed_sequences.fasta"
    if not passed_fasta.exists():
        print("\n✗ No sequences passed QC. Pipeline stopped.")
        sys.exit(1)

    # Step 2: Species Identification
    id_script = Path(__file__).parent / "04_identification" / "identify_species.py"
    id_success = run_command(
        ["python", str(id_script), str(passed_fasta), str(run_dir)],
        "Species Identification"
    )

    if not id_success:
        print("\n⚠ Pipeline completed with warnings (identification failed)")

    # Final summary
    print(f"\n{'='*60}")
    print("Pipeline Complete!")
    print(f"{'='*60}")
    print(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"\nResults saved to: {run_dir}")
    print("\nGenerated files:")
    print(f"  - qc_report.html - Quality control summary")
    print(f"  - passed_sequences.fasta - Sequences that passed QC")
    if id_success:
        print(f"  - identification_report.html - Species identification results")

    print("\nNext steps:")
    print(f"  1. Open {run_dir}/qc_report.html in your browser")
    if id_success:
        print(f"  2. Review species identifications in identification_report.html")

if __name__ == "__main__":
    main()
