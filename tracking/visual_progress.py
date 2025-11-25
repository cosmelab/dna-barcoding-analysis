#!/usr/bin/env python3
"""
Visual Progress Display for DNA Barcoding Analysis Repository
Creates a clear visual representation of project progress
Uses tracking v3.0 (master_state.json format)
"""

import json
from pathlib import Path
from datetime import datetime

# Try to import rich for pretty output
try:
    from rich.console import Console
    from rich.panel import Panel
    from rich.table import Table
    from rich.text import Text
    RICH_AVAILABLE = True
    console = Console()
except ImportError:
    RICH_AVAILABLE = False
    console = None

PROJECT_ROOT = Path(__file__).parent.parent
TRACKING_DIR = PROJECT_ROOT / "tracking"
STATE_FILE = TRACKING_DIR / "master_state.json"

def load_state():
    if not STATE_FILE.exists():
        return {
            "cycle": 0,
            "last_updated": "Not initialized",
            "current_session": {"focus": "Tracking system not yet initialized"},
            "active_tasks": [],
            "completed_this_cycle": []
        }
    with open(STATE_FILE, 'r') as f:
        return json.load(f)

def create_visual():
    state = load_state()

    if RICH_AVAILABLE:
        create_visual_rich(state)
    else:
        create_visual_plain(state)


def create_visual_rich(state):
    """Create visual progress with rich formatting"""
    # Header
    console.print(Panel.fit(
        "[bold cyan]DNA BARCODING ANALYSIS[/bold cyan]\n[dim]ENTM201L GitHub Classroom Repository[/dim]",
        border_style="cyan"
    ))

    # Info table
    info_table = Table(show_header=False, box=None, padding=(0, 2))
    info_table.add_column("Key", style="dim")
    info_table.add_column("Value")
    info_table.add_row("Generated", datetime.now().strftime('%Y-%m-%d %H:%M'))
    info_table.add_row("Cycle", str(state.get('cycle', 0)))
    info_table.add_row("Last Updated", state.get('last_updated', 'Not initialized'))
    console.print(info_table)

    # Repository Purpose
    console.print("\n[bold cyan]REPOSITORY PURPOSE[/bold cyan]")
    console.print("[dim]─" * 60 + "[/dim]")
    console.print("Teaching DNA barcoding bioinformatics for ENTM201L students")
    console.print("Two pathways: GUI (beginner-friendly) and CLI (reproducible)")

    # Analysis Pathways
    console.print("\n[bold cyan]ANALYSIS PATHWAYS[/bold cyan]")
    console.print("[dim]─" * 60 + "[/dim]")

    pathway_table = Table(show_header=True, header_style="bold")
    pathway_table.add_column("GUI Pathway", style="green")
    pathway_table.add_column("CLI Pathway", style="blue")
    pathway_table.add_row("Tracy web → QC", "Python scripts → QC")
    pathway_table.add_row("UGENE/SeaView → Align", "MAFFT → Align")
    pathway_table.add_row("Built-in trees", "IQ-TREE2 → Trees")
    pathway_table.add_row("[dim]No Docker needed[/dim]", "[dim]Docker container[/dim]")
    console.print(pathway_table)

    # Module Status
    console.print("\n[bold cyan]MODULE STATUS[/bold cyan]")
    console.print("[dim]─" * 60 + "[/dim]")

    modules = [
        ("modules/01_quality_control", "Quality Control"),
        ("modules/02_consensus", "Consensus Generation"),
        ("modules/03_alignment", "Sequence Alignment"),
        ("modules/04_phylogeny", "Phylogenetic Trees"),
        ("modules/05_identification", "Species ID (BLAST)"),
    ]

    for module_path, module_name in modules:
        module_dir = PROJECT_ROOT / module_path
        if module_dir.exists():
            console.print(f"  [green]✓[/green] {module_name}")
        else:
            console.print(f"  [dim]○[/dim] {module_name}")

    # Active Tasks
    active_tasks = state.get('active_tasks', [])
    console.print("\n[bold cyan]ACTIVE TASKS[/bold cyan]")
    console.print("[dim]─" * 60 + "[/dim]")

    if active_tasks:
        for task in active_tasks:
            status_color = {'pending': 'yellow', 'in_progress': 'blue', 'blocked': 'red'}.get(task.get('status'), 'white')
            priority = " [red][HIGH][/red]" if task.get('priority') == 'HIGH' else ""
            console.print(f"  [{status_color}]●[/{status_color}] {task['task']}{priority}")
    else:
        console.print("  [green]✓[/green] No active tasks - ready for students!")

    # Development Status
    console.print("\n[bold cyan]STATUS[/bold cyan]")
    console.print("[dim]─" * 60 + "[/dim]")
    console.print("  [green]✓[/green] GitHub Classroom: Ready")
    console.print("  [green]✓[/green] Docker container: Available")
    console.print("  [green]✓[/green] All modules: Complete")

    # Commands
    console.print("\n[bold cyan]COMMANDS[/bold cyan]")
    console.print("[dim]─" * 60 + "[/dim]")
    console.print("  [dim]python3 tracking/monitor.py status[/dim]")
    console.print("  [dim]python3 tracking/visual_progress.py[/dim]")
    console.print()

    # Also save plain text version
    save_plain_output(state)


def save_plain_output(state):
    """Save plain text version to file"""
    output = create_plain_output(state)
    output_file = TRACKING_DIR / 'visual_progress.txt'
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        f.write('\n'.join(output))
    console.print(f"[dim]Saved to: {output_file}[/dim]")


def create_plain_output(state):
    """Generate plain text output lines"""
    output = []
    output.append("="*80)
    output.append("DNA BARCODING ANALYSIS - ENTM201L GITHUB CLASSROOM REPO")
    output.append("="*80)
    output.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    output.append(f"Cycle: {state.get('cycle', 0)}")
    output.append(f"Last Updated: {state.get('last_updated', 'Not initialized')}")
    output.append("="*80)
    output.append("\nMODULE STATUS:")
    output.append("-" * 80)

    modules = [
        ("modules/01_quality_control", "Quality Control"),
        ("modules/02_consensus", "Consensus Generation"),
        ("modules/03_alignment", "Sequence Alignment"),
        ("modules/04_phylogeny", "Phylogenetic Trees"),
        ("modules/05_identification", "Species ID (BLAST)"),
    ]

    for module_path, module_name in modules:
        module_dir = PROJECT_ROOT / module_path
        status = "✓" if module_dir.exists() else "○"
        output.append(f"  {status} {module_name}")

    output.append("\n" + "="*80)
    return output


def create_visual_plain(state):
    """Create visual progress with plain text (fallback)"""
    output = create_plain_output(state)

    # Write to file
    output_file = TRACKING_DIR / 'visual_progress.txt'
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        f.write('\n'.join(output))

    # Print to console
    print('\n'.join(output))
    print(f"\n✅ Visual progress saved to: {output_file}")

if __name__ == "__main__":
    create_visual()
