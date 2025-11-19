#!/usr/bin/env python3
"""
Visual Progress Display for DNA Barcoding Analysis Repository
Creates a clear visual representation of project progress
Uses tracking v3.0 (master_state.json format)
"""

import json
from pathlib import Path
from datetime import datetime

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

    output = []

    # Header
    output.append("="*80)
    output.append("DNA BARCODING ANALYSIS - ENTM201L GITHUB CLASSROOM REPO")
    output.append("="*80)
    output.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    output.append(f"Cycle: {state.get('cycle', 0)}")
    output.append(f"Last Updated: {state.get('last_updated', 'Not initialized')}")
    output.append("="*80)

    # Repository Purpose
    output.append("\n📚 REPOSITORY PURPOSE:")
    output.append("-" * 80)
    output.append("Teaching DNA barcoding bioinformatics for ENTM201L students")
    output.append("Two pathways: GUI (beginner-friendly) and CLI (reproducible science)")
    output.append("GitHub Classroom ready - students get individual repositories")

    # Current Session Focus
    session = state.get('current_session', {})
    if session:
        output.append("\n🎯 CURRENT FOCUS:")
        output.append("-" * 80)
        output.append(f"Session: {session.get('focus', 'Unknown')}")
        output.append(f"Environment: {session.get('environment', 'Unknown')}")
        output.append(f"Status: {session.get('user_status', 'Unknown')}")

    # Analysis Pathways
    output.append("\n🔬 ANALYSIS PATHWAYS:")
    output.append("-" * 80)
    output.append("GUI PATHWAY (Easiest - No command-line required):")
    output.append("  • Tracy web interface → chromatogram QC")
    output.append("  • UGENE or SeaView → alignment + phylogenetic trees")
    output.append("  • No Docker required")
    output.append("")
    output.append("CLI PATHWAY (Reproducible - Learn bioinformatics):")
    output.append("  • Tracy CLI → chromatogram QC")
    output.append("  • MAFFT → sequence alignment")
    output.append("  • IQ-TREE2 → phylogenetic trees")
    output.append("  • ggtree/FigTree → tree visualization")
    output.append("  • All in Docker container")

    # Active Tasks
    output.append("\n📋 ACTIVE TASKS:")
    output.append("-" * 80)

    active_tasks = state.get('active_tasks', [])
    if active_tasks:
        for task in active_tasks:
            status_icon = {
                'pending': '⏳',
                'in_progress': '🔄',
                'blocked': '🚫'
            }.get(task.get('status', 'pending'), '❓')

            priority_badge = ''
            if task.get('priority') == 'HIGH':
                priority_badge = ' [HIGH PRIORITY]'

            output.append(f"{status_icon} {task['task']}{priority_badge}")
            output.append(f"   ID: {task['id']}")
            if task.get('description'):
                output.append(f"   {task['description']}")
            output.append("")
    else:
        output.append("✅ No active tasks - ready for students!")

    # Module Status
    output.append("\n📖 MODULE STATUS:")
    output.append("-" * 80)
    modules = [
        "00_introduction - Course overview and setup",
        "01_linux_basics - Command-line fundamentals",
        "02_python_basics - Python for bioinformatics",
        "03_r_basics - R and data visualization",
        "04_data - Sample sequencing data",
        "05_quality_control - Sanger chromatogram QC",
        "06_alignment - Multiple sequence alignment",
        "07_phylogeny - Phylogenetic tree building",
        "08_identification - Species identification"
    ]

    for module in modules:
        module_dir = PROJECT_ROOT / module.split(' - ')[0]
        if module_dir.exists():
            output.append(f"  ✅ {module}")
        else:
            output.append(f"  ⏳ {module}")

    # Completed This Cycle
    completed = state.get('completed_this_cycle', [])
    if completed:
        output.append(f"\n✅ COMPLETED THIS CYCLE ({len(completed)} tasks):")
        output.append("-" * 80)
        for task in completed[-5:]:  # Show last 5
            output.append(f"  • {task.get('task', 'Unknown')}")
            if task.get('completed'):
                completed_date = task['completed'].split('T')[0]
                output.append(f"    Completed: {completed_date}")

    # Development Status
    output.append("\n🚧 DEVELOPMENT STATUS:")
    output.append("-" * 80)
    output.append("⚠️  Repository is in active development")
    output.append("⚠️  Students: Do NOT run scripts yet")
    output.append("✅ Expected ready: Week 8 (November 20, 2025)")
    output.append("✅ GitHub Classroom integration: Ready")
    output.append("✅ Docker container: Available")

    # Key Resources
    output.append("\n🔗 KEY RESOURCES:")
    output.append("-" * 80)
    output.append("  • Setup guide: ../entm201l-fall2025/setup/docker-installation.html")
    output.append("  • Course website: http://138.23.14.176:8001/")
    output.append("  • Container: cosmelab/entm201l:latest")
    output.append("  • Tracy web QC: https://www.gear-genomics.com")

    # Footer
    output.append("\n" + "="*80)
    output.append("COMMANDS:")
    output.append("  • View status: python3 tracking/monitor.py status")
    output.append("  • Update progress: python3 tracking/visual_progress.py")
    output.append("="*80)

    # Write to file
    output_file = TRACKING_DIR / 'visual_progress.txt'
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        f.write('\n'.join(output))

    # Also print to console
    print('\n'.join(output))
    print(f"\n✅ Visual progress saved to: {output_file}")

if __name__ == "__main__":
    create_visual()
