#!/usr/bin/env python3
"""
DNA Barcoding Master Pipeline
Runs complete analysis: Assembly -> QC -> Alignment -> Phylogeny -> Identification
"""

import os
import sys
import subprocess
import argparse
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
    except FileNotFoundError:
        print(f"✗ {step_name} failed - script not found")
        return False

def main():
    """Run complete DNA barcoding pipeline"""
    parser = argparse.ArgumentParser(
        description="DNA Barcoding Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Full pipeline (recommended)
  python master_pipeline.py data/my_sequences/

  # Skip assembly (if you already have merged sequences)
  python master_pipeline.py data/my_sequences/ --skip-assembly

  # Skip phylogeny (faster, ID only)
  python master_pipeline.py data/my_sequences/ --skip-phylogeny

  # Reference sequences for alignment
  python master_pipeline.py data/my_sequences/ --reference data/reference/socal_mosquitoes.fasta
        """)

    parser.add_argument("input_dir", help="Directory containing .ab1 chromatogram files")
    parser.add_argument("output_dir", nargs="?", default="results", help="Output directory (default: results/)")
    parser.add_argument("--skip-assembly", action="store_true", help="Skip F/R sequence assembly step")
    parser.add_argument("--skip-alignment", action="store_true", help="Skip alignment step")
    parser.add_argument("--skip-phylogeny", action="store_true", help="Skip phylogeny step")
    parser.add_argument("--reference", help="Reference sequences FASTA for alignment")

    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

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
    print(f"\nPipeline steps:")
    if not args.skip_assembly:
        print("  1. Assembly (F/R merging)")
    print(f"  {'2' if not args.skip_assembly else '1'}. Quality Control")
    if not args.skip_alignment:
        print(f"  {'3' if not args.skip_assembly else '2'}. Alignment")
    if not args.skip_phylogeny:
        print(f"  {'4' if not args.skip_assembly else '3'}. Phylogeny")
    print(f"  {'5' if not args.skip_assembly else '4'}. Species Identification")

    # Module base path
    module_dir = Path(__file__).parent

    # Working FASTA file - will be updated through pipeline
    working_fasta = None

    # Step 0: Assembly (F/R merging) - OPTIONAL
    if not args.skip_assembly:
        assembly_script = module_dir / "00_assembly" / "merge_forward_reverse.py"
        assembly_dir = run_dir / "00_assembly"
        assembly_dir.mkdir(exist_ok=True)

        assembly_success = run_command(
            ["python", str(assembly_script), str(input_dir), str(assembly_dir)],
            "F/R Sequence Assembly"
        )

        if assembly_success:
            working_fasta = assembly_dir / "consensus_sequences.fasta"
        else:
            print("\n⚠ Assembly failed - continuing with individual sequences")

    # Step 1: Quality Control
    qc_script = module_dir / "01_quality_control" / "qc_chromatograms.py"
    qc_dir = run_dir / "01_qc"
    qc_dir.mkdir(exist_ok=True)

    qc_success = run_command(
        ["python", str(qc_script), str(input_dir), str(qc_dir)],
        "Quality Control"
    )

    if not qc_success:
        print("\n✗ Pipeline failed at Quality Control step")
        sys.exit(1)

    # Check if we have passed sequences
    passed_fasta = qc_dir / "passed_sequences.fasta"
    if not passed_fasta.exists():
        print("\n✗ No sequences passed QC. Pipeline stopped.")
        sys.exit(1)

    # Update working FASTA to QC-passed sequences
    working_fasta = passed_fasta

    # Step 2: Alignment (OPTIONAL)
    aligned_fasta = None
    if not args.skip_alignment:
        alignment_script = module_dir / "02_alignment" / "align_sequences.py"
        alignment_dir = run_dir / "02_alignment"
        alignment_dir.mkdir(exist_ok=True)

        # Combine with reference if provided
        alignment_input = working_fasta
        if args.reference:
            ref_path = Path(args.reference)
            if ref_path.exists():
                print(f"\nCombining with reference sequences from {ref_path}")
                combined = alignment_dir / "input_with_reference.fasta"
                # Combine sequences
                with open(combined, 'w') as out:
                    with open(working_fasta) as f:
                        out.write(f.read())
                    with open(ref_path) as f:
                        out.write(f.read())
                alignment_input = combined

        alignment_success = run_command(
            ["python", str(alignment_script), str(alignment_input), str(alignment_dir)],
            "Sequence Alignment"
        )

        if alignment_success:
            aligned_fasta = alignment_dir / "aligned_sequences.fasta"
            working_fasta = aligned_fasta
        else:
            print("\n⚠ Alignment failed - continuing with unaligned sequences")

    # Step 3: Phylogeny (OPTIONAL)
    if not args.skip_phylogeny and aligned_fasta:
        phylogeny_script = module_dir / "03_phylogeny" / "build_tree.py"
        phylogeny_dir = run_dir / "03_phylogeny"
        phylogeny_dir.mkdir(exist_ok=True)

        phylogeny_success = run_command(
            ["python", str(phylogeny_script), str(aligned_fasta), str(phylogeny_dir)],
            "Phylogenetic Tree"
        )

        if not phylogeny_success:
            print("\n⚠ Phylogeny failed - continuing to species ID")

    # Step 4: Species Identification
    id_script = module_dir / "04_identification" / "identify_species.py"
    id_dir = run_dir / "04_identification"
    id_dir.mkdir(exist_ok=True)

    id_success = run_command(
        ["python", str(id_script), str(passed_fasta), str(id_dir)],
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
    print(f"  01_qc/qc_report.html - Quality control summary")
    print(f"  01_qc/passed_sequences.fasta - Sequences that passed QC")
    if not args.skip_assembly:
        print(f"  00_assembly/consensus_sequences.fasta - Merged F/R sequences")
    if not args.skip_alignment and aligned_fasta:
        print(f"  02_alignment/aligned_sequences.fasta - Multiple sequence alignment")
    if not args.skip_phylogeny and aligned_fasta:
        print(f"  03_phylogeny/tree.png - Phylogenetic tree visualization")
        print(f"  03_phylogeny/tree.treefile - Tree in Newick format (open with FigTree)")
    if id_success:
        print(f"  04_identification/identification_report.html - Species identification results")

    print("\nNext steps:")
    print(f"  1. Open {run_dir}/01_qc/qc_report.html in your browser")
    if not args.skip_phylogeny:
        print(f"  2. View phylogenetic tree: {run_dir}/03_phylogeny/tree.png")
    if id_success:
        print(f"  3. Review species identifications: {run_dir}/04_identification/identification_report.html")

if __name__ == "__main__":
    main()
