#!/usr/bin/env python3
"""
Student Progress Tracker for DNA Barcoding Analysis
Automatically detects which pipeline steps have been completed.

Usage:
    python3 check_progress.py          # Check progress (works locally or in container)

Inside Docker container:
    docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
      cosmelab/dna-barcoding-analysis:latest python3 check_progress.py
"""

from pathlib import Path
import os

# Try to import rich for pretty output
try:
    from rich.console import Console
    from rich.table import Table
    from rich.panel import Panel
    from rich.progress import Progress, BarColumn, TextColumn
    RICH_AVAILABLE = True
    # force_terminal=True ensures colors work even when piped through tee
    console = Console(force_terminal=True)
except ImportError:
    RICH_AVAILABLE = False
    console = None

# Detect if running in container or locally
if Path('/workspace').exists() and Path('/workspace/modules').exists():
    PROJECT_ROOT = Path('/workspace')
else:
    PROJECT_ROOT = Path(__file__).parent


def check_step_completion(results_dir: Path) -> dict:
    """Check which pipeline steps are complete based on output files."""
    steps = {
        'qc': {
            'name': '1. Quality Control',
            'files': ['qc_report.html', 'passed_sequences.fasta'],
            'dir': results_dir / '01_qc',
            'complete': False
        },
        'consensus': {
            'name': '2. Consensus Sequences',
            'files': ['consensus_sequences.fasta', 'consensus_report.html'],
            'dir': results_dir / '02_consensus',
            'complete': False
        },
        'alignment': {
            'name': '3. Sequence Alignment',
            'files': ['aligned_sequences.fasta', 'alignment_report.html'],
            'dir': results_dir / '03_alignment',
            'complete': False
        },
        'phylogeny': {
            'name': '4. Phylogenetic Tree',
            'files': ['tree.treefile', 'phylogeny_report.html'],
            'dir': results_dir / '04_phylogeny',
            'complete': False
        },
        'blast': {
            'name': '5. Species ID (BLAST)',
            'files': ['identification_report.html'],
            'dir': results_dir / '05_blast',
            'complete': False
        }
    }

    for step_id, step in steps.items():
        step_dir = step['dir']
        if step_dir.exists():
            # Check if at least one required file exists
            files_found = [f for f in step['files'] if (step_dir / f).exists()]
            step['complete'] = len(files_found) > 0
            step['files_found'] = files_found
        else:
            step['files_found'] = []

    return steps


def calculate_progress(steps: dict) -> tuple:
    """Calculate completion percentage."""
    completed = sum(1 for s in steps.values() if s['complete'])
    total = len(steps)
    return completed, total, int((completed / total) * 100)


def display_progress_rich():
    """Display progress using rich formatting."""
    tutorial_dir = PROJECT_ROOT / 'results' / 'tutorial'
    my_analysis_dir = PROJECT_ROOT / 'results' / 'my_analysis'

    tutorial_steps = check_step_completion(tutorial_dir)
    my_analysis_steps = check_step_completion(my_analysis_dir)

    tut_done, tut_total, tut_pct = calculate_progress(tutorial_steps)
    my_done, my_total, my_pct = calculate_progress(my_analysis_steps)

    # Header
    console.print()
    console.print(Panel.fit(
        "[bold cyan]DNA Barcoding Analysis - Progress Tracker[/bold cyan]",
        border_style="cyan"
    ))

    # Progress bars
    console.print()
    console.print(f"[bold]Tutorial Progress:[/bold]    ", end="")
    bar_tut = "[green]" + ("█" * tut_done) + "[/green][dim]" + ("░" * (tut_total - tut_done)) + "[/dim]"
    console.print(f"{bar_tut} {tut_done}/{tut_total} ({tut_pct}%)")

    console.print(f"[bold]Your Analysis:[/bold]       ", end="")
    bar_my = "[blue]" + ("█" * my_done) + "[/blue][dim]" + ("░" * (my_total - my_done)) + "[/dim]"
    console.print(f"{bar_my} {my_done}/{my_total} ({my_pct}%)")

    # Detailed table
    console.print()
    table = Table(title="Step-by-Step Status", show_header=True, header_style="bold")
    table.add_column("Pipeline Step", style="white", width=25)
    table.add_column("Tutorial", justify="center", width=10)
    table.add_column("Your Data", justify="center", width=10)

    step_keys = ['qc', 'consensus', 'alignment', 'phylogeny', 'blast']
    for key in step_keys:
        tut_status = "[green]✓[/green]" if tutorial_steps[key]['complete'] else "[dim]○[/dim]"
        my_status = "[green]✓[/green]" if my_analysis_steps[key]['complete'] else "[dim]○[/dim]"
        table.add_row(tutorial_steps[key]['name'], tut_status, my_status)

    console.print(table)

    # Next steps hint
    console.print()
    if tut_done < tut_total:
        next_step = next(s for s in tutorial_steps.values() if not s['complete'])
        console.print(f"[yellow]Next tutorial step:[/yellow] {next_step['name']}")
        console.print("[dim]Run: ./tutorial.sh (or see README for manual commands)[/dim]")
    elif my_done < my_total:
        next_step = next(s for s in my_analysis_steps.values() if not s['complete'])
        console.print(f"[yellow]Next step for your data:[/yellow] {next_step['name']}")
        console.print("[dim]Run: ./run-analysis.sh (or see README for manual commands)[/dim]")
    else:
        console.print("[bold green]All steps complete! Review your results and fill out assignment.md[/bold green]")

    console.print()


def display_progress_plain():
    """Display progress using plain text (fallback)."""
    tutorial_dir = PROJECT_ROOT / 'results' / 'tutorial'
    my_analysis_dir = PROJECT_ROOT / 'results' / 'my_analysis'

    tutorial_steps = check_step_completion(tutorial_dir)
    my_analysis_steps = check_step_completion(my_analysis_dir)

    tut_done, tut_total, tut_pct = calculate_progress(tutorial_steps)
    my_done, my_total, my_pct = calculate_progress(my_analysis_steps)

    print()
    print("=" * 60)
    print("DNA BARCODING ANALYSIS - PROGRESS TRACKER")
    print("=" * 60)

    # Progress bars
    bar_tut = ("█" * tut_done) + ("░" * (tut_total - tut_done))
    bar_my = ("█" * my_done) + ("░" * (my_total - my_done))

    print(f"\nTutorial Progress:    {bar_tut} {tut_done}/{tut_total} ({tut_pct}%)")
    print(f"Your Analysis:        {bar_my} {my_done}/{my_total} ({my_pct}%)")

    print("\n" + "-" * 60)
    print(f"{'Pipeline Step':<25} {'Tutorial':^10} {'Your Data':^10}")
    print("-" * 60)

    step_keys = ['qc', 'consensus', 'alignment', 'phylogeny', 'blast']
    for key in step_keys:
        tut_status = "✓" if tutorial_steps[key]['complete'] else "○"
        my_status = "✓" if my_analysis_steps[key]['complete'] else "○"
        print(f"{tutorial_steps[key]['name']:<25} {tut_status:^10} {my_status:^10}")

    print("-" * 60)

    # Next steps
    if tut_done < tut_total:
        next_step = next(s for s in tutorial_steps.values() if not s['complete'])
        print(f"\nNext tutorial step: {next_step['name']}")
    elif my_done < my_total:
        next_step = next(s for s in my_analysis_steps.values() if not s['complete'])
        print(f"\nNext step for your data: {next_step['name']}")
    else:
        print("\nAll steps complete! Review your results and fill out assignment.md")

    print()


def main():
    """Main entry point."""
    if RICH_AVAILABLE:
        display_progress_rich()
    else:
        display_progress_plain()


if __name__ == "__main__":
    main()
