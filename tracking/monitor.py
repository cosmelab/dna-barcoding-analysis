#!/usr/bin/env python3
"""
Tracking Monitor v3.0 - Self-Maintaining Tracking System

This script monitors and enforces the tracking system rules:
- Checks if cycle should be archived (10 tasks or 24 hours)
- Consolidates scattered tracking data
- Enforces file naming and count rules
- Provides automatic archiving functionality
"""

import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple, Any

# Try to import rich for pretty output, fall back to plain text
try:
    from rich.console import Console
    from rich.panel import Panel
    from rich.table import Table
    from rich.progress import Progress, BarColumn, TextColumn
    from rich import print as rprint
    RICH_AVAILABLE = True
    console = Console()
except ImportError:
    RICH_AVAILABLE = False
    console = None


class TrackingMonitor:
    """Monitor and enforce tracking system v3.0 rules"""

    TRACKING_DIR = Path(__file__).parent
    MASTER_STATE = TRACKING_DIR / "master_state.json"
    ARCHIVE_INDEX = TRACKING_DIR / "archive_index.json"
    ARCHIVE_DIR = TRACKING_DIR / "archive"

    # Allowed files in tracking/ directory (ONLY these!)
    ALLOWED_FILES = {
        "master_state.json",
        "archive_index.json",
        "CHANGELOG.md",
        "STYLE_GUIDE.md",  # Special case: ONLY allowed ALL CAPS file
        "style_guide.md",  # Lowercase version also allowed
        "assistant_rules.md",
        "tracking_v3_design.md",
        "monitor.py",
        "__pycache__"  # Python cache directory
    }

    def __init__(self):
        """Initialize monitor by loading current state"""
        self.state = self._load_json(self.MASTER_STATE)
        self.index = self._load_json(self.ARCHIVE_INDEX)

    def _load_json(self, filepath: Path) -> Dict[str, Any]:
        """Load JSON file with error handling"""
        if not filepath.exists():
            print(f"ERROR: {filepath} not found!")
            return {}

        try:
            with open(filepath, 'r') as f:
                return json.load(f)
        except json.JSONDecodeError as e:
            print(f"ERROR: Invalid JSON in {filepath}: {e}")
            return {}

    def _save_json(self, filepath: Path, data: Dict[str, Any]) -> None:
        """Save JSON file with pretty formatting"""
        with open(filepath, 'w') as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        print(f"✓ Updated {filepath.name}")

    def check_cycle(self) -> Tuple[str, str]:
        """
        Check if current cycle should be archived

        Returns:
            Tuple of (action, reason) where action is 'archive' or 'continue'
        """
        # Get completed tasks count
        completed = self.state.get("completed_this_cycle", [])
        task_count = len(completed)

        # Calculate time elapsed
        cycle_started = self.state.get("cycle_started", "")
        if cycle_started:
            start_time = datetime.fromisoformat(cycle_started.replace('Z', '+00:00'))
            now = datetime.now(timezone.utc)
            elapsed_hours = (now - start_time).total_seconds() / 3600
        else:
            elapsed_hours = 0

        # Check trigger conditions
        if task_count >= 10:
            return ("archive", f"10 tasks completed ({task_count})")

        if elapsed_hours >= 24:
            return ("archive", f"24 hours elapsed ({elapsed_hours:.1f}h)")

        return ("continue", None)

    def enforce_rules(self) -> List[Dict[str, Any]]:
        """
        Check for tracking system violations

        Returns:
            List of violation dictionaries
        """
        violations = []

        # Check for extra tracking files
        tracking_files = [
            f.name for f in self.TRACKING_DIR.iterdir()
            if f.is_file() and f.suffix in ['.json', '.md']
        ]

        extra_files = [
            f for f in tracking_files
            if f not in self.ALLOWED_FILES
        ]

        if extra_files:
            violations.append({
                "type": "extra_files",
                "files": extra_files,
                "action": "move_to_archive",
                "severity": "medium"
            })

        # Check for ALL CAPS filenames (except STYLE_GUIDE.md)
        all_files = [f.name for f in self.TRACKING_DIR.iterdir() if f.is_file()]
        caps_files = [
            f for f in all_files
            if f.isupper() and f not in ["STYLE_GUIDE.md", "CHANGELOG.md"]
        ]

        if caps_files:
            violations.append({
                "type": "caps_filenames",
                "files": caps_files,
                "action": "rename_to_lowercase",
                "severity": "high"
            })

        # Check for directories that shouldn't exist
        unwanted_dirs = [
            d.name for d in self.TRACKING_DIR.iterdir()
            if d.is_dir() and d.name not in ["archive", "__pycache__", "citation_reports"]
        ]

        if unwanted_dirs:
            violations.append({
                "type": "unwanted_directories",
                "directories": unwanted_dirs,
                "action": "review_and_archive",
                "severity": "low"
            })

        return violations

    def consolidate(self) -> Dict[str, Any]:
        """
        Consolidate any scattered tracking data into master_state.json

        Returns:
            Dictionary with consolidation results
        """
        results = {
            "files_checked": [],
            "data_merged": [],
            "warnings": []
        }

        # Check for agent_reports.json (shouldn't exist in v3.0)
        agent_reports = self.TRACKING_DIR / "agent_reports.json"
        if agent_reports.exists():
            results["warnings"].append(
                "agent_reports.json exists - should be in archive"
            )
            results["files_checked"].append("agent_reports.json")

        # Check for master_tasks.json (v2.0 format, should be migrated)
        master_tasks = self.TRACKING_DIR / "master_tasks.json"
        if master_tasks.exists():
            results["warnings"].append(
                "master_tasks.json exists - v2.0 format should be archived"
            )
            results["files_checked"].append("master_tasks.json")

        return results

    def get_status(self) -> Dict[str, Any]:
        """Get current tracking system status"""
        action, reason = self.check_cycle()
        violations = self.enforce_rules()
        consolidation = self.consolidate()

        # Calculate metrics
        active_tasks = len(self.state.get("active_tasks", []))
        completed_tasks = len(self.state.get("completed_this_cycle", []))
        cycle = self.state.get("cycle", 0)

        return {
            "tracking_version": "3.0",
            "cycle": cycle,
            "active_tasks": active_tasks,
            "completed_tasks": completed_tasks,
            "archive_action": action,
            "archive_reason": reason,
            "violations": violations,
            "consolidation": consolidation,
            "health": "good" if not violations else "warnings"
        }

    def print_status(self) -> None:
        """Print human-readable status report"""
        status = self.get_status()

        if RICH_AVAILABLE:
            self._print_status_rich(status)
        else:
            self._print_status_plain(status)

    def _print_status_rich(self, status: Dict[str, Any]) -> None:
        """Print status using rich formatting"""
        # Header panel
        console.print(Panel.fit(
            "[bold cyan]TRACKING SYSTEM v3.0[/bold cyan]",
            border_style="cyan"
        ))

        # Status table
        table = Table(show_header=False, box=None, padding=(0, 2))
        table.add_column("Key", style="dim")
        table.add_column("Value")

        table.add_row("Cycle", f"[bold]{status['cycle']}[/bold]")
        table.add_row("Active tasks", str(status['active_tasks']))

        # Progress bar for completed tasks
        completed = status['completed_tasks']
        progress_bar = "█" * completed + "░" * (10 - completed)
        table.add_row("Completed", f"[green]{progress_bar}[/green] {completed}/10")

        # Archive status with color
        archive_color = "yellow" if status['archive_action'] == "archive" else "green"
        table.add_row("Archive", f"[{archive_color}]{status['archive_action'].upper()}[/{archive_color}]")

        if status['archive_reason']:
            table.add_row("", f"[dim]→ {status['archive_reason']}[/dim]")

        # Health status
        health_color = "green" if status['health'] == "good" else "yellow"
        table.add_row("Health", f"[{health_color}]{status['health'].upper()}[/{health_color}]")

        console.print(table)

        # Violations
        if status['violations']:
            console.print(f"\n[yellow]⚠ VIOLATIONS ({len(status['violations'])})[/yellow]")
            for v in status['violations']:
                items = v.get('files', v.get('directories', []))
                console.print(f"  [red]•[/red] {v['type']}: {len(items)} items")
                console.print(f"    [dim]Action: {v['action']}[/dim]")

        # Consolidation warnings
        if status['consolidation']['warnings']:
            console.print(f"\n[yellow]⚠ CONSOLIDATION WARNINGS[/yellow]")
            for w in status['consolidation']['warnings']:
                console.print(f"  [yellow]•[/yellow] {w}")

        console.print()

    def _print_status_plain(self, status: Dict[str, Any]) -> None:
        """Print status using plain text (fallback)"""
        print("\n" + "="*60)
        print("TRACKING SYSTEM v3.0 STATUS")
        print("="*60)
        print(f"Cycle: {status['cycle']}")
        print(f"Active tasks: {status['active_tasks']}")
        print(f"Completed this cycle: {status['completed_tasks']}/10")
        print(f"Archive status: {status['archive_action'].upper()}")
        if status['archive_reason']:
            print(f"  → Reason: {status['archive_reason']}")

        print(f"\nHealth: {status['health'].upper()}")

        if status['violations']:
            print(f"\n⚠️  VIOLATIONS DETECTED ({len(status['violations'])})")
            for v in status['violations']:
                print(f"  • {v['type']}: {len(v.get('files', v.get('directories', [])))} items")
                print(f"    Action: {v['action']}")

        if status['consolidation']['warnings']:
            print(f"\n⚠️  CONSOLIDATION WARNINGS")
            for w in status['consolidation']['warnings']:
                print(f"  • {w}")

        print("="*60 + "\n")

    def archive_cycle(self, summary: str = None) -> bool:
        """
        Archive current cycle and start new one

        Args:
            summary: Human-readable summary slug (auto-generated if not provided)

        Returns:
            True if successful, False otherwise
        """
        # Generate timestamp and summary
        now = datetime.now(timezone.utc)
        timestamp = now.strftime("%Y-%m-%d_%H-%M")

        if not summary:
            # Auto-generate summary from milestones or focus
            milestones = self.state.get("milestones_this_cycle", [])
            if milestones:
                summary = milestones[-1]["name"].lower().replace(" ", "-")[:30]
            else:
                focus = self.state.get("current_session", {}).get("focus", "work")
                summary = focus.lower().replace(" ", "-")[:30]

        # Create archive filename
        archive_file = self.ARCHIVE_DIR / f"{timestamp}_{summary}.json"

        # Save current state to archive
        try:
            self._save_json(archive_file, self.state)
            print(f"✓ Archived Cycle {self.state['cycle']} to {archive_file.name}")
        except Exception as e:
            print(f"ERROR: Failed to create archive: {e}")
            return False

        # Update archive index
        cycle_num = self.state.get("cycle", 0)
        cycle_started = self.state.get("cycle_started", "")

        archive_entry = {
            "cycle": cycle_num,
            "file": f"archive/{archive_file.name}",
            "started": cycle_started,
            "ended": now.isoformat(),
            "duration": self._calculate_duration(cycle_started, now.isoformat()),
            "tasks_completed": len(self.state.get("completed_this_cycle", [])),
            "milestones": [m["name"] for m in self.state.get("milestones_this_cycle", [])],
            "summary": summary.replace("-", " ").title()
        }

        self.index["archives"].append(archive_entry)
        self.index["total_cycles"] = len(self.index["archives"])
        self.index["total_tasks_completed"] = sum(
            a["tasks_completed"] for a in self.index["archives"]
        )

        self._save_json(self.ARCHIVE_INDEX, self.index)

        # Reset master_state.json for new cycle
        new_state = {
            "tracking_version": "3.0",
            "cycle": cycle_num + 1,
            "cycle_started": now.isoformat(),
            "last_updated": now.isoformat(),
            "archive_trigger": "task_count_10_or_24_hours",
            "current_session": self.state.get("current_session", {}),
            "active_tasks": [],
            "completed_this_cycle": [],
            "milestones_this_cycle": [],
            "files_modified_this_cycle": [],
            "next_archive_trigger": {
                "reason": "10 tasks completed OR 24 hours elapsed",
                "current_task_count": 0,
                "current_time_elapsed": "0h 0min",
                "will_archive_when": "10 tasks OR 24 hours"
            },
            "context": self.state.get("context", {})
        }

        self._save_json(self.MASTER_STATE, new_state)

        print(f"✓ Started Cycle {cycle_num + 1}")
        print(f"\n🎉 Cycle {cycle_num} complete!")

        return True

    def _calculate_duration(self, start: str, end: str) -> str:
        """Calculate human-readable duration between timestamps"""
        try:
            start_dt = datetime.fromisoformat(start.replace('Z', '+00:00'))
            end_dt = datetime.fromisoformat(end.replace('Z', '+00:00'))
            delta = end_dt - start_dt

            hours = int(delta.total_seconds() / 3600)
            minutes = int((delta.total_seconds() % 3600) / 60)

            if hours > 0:
                return f"{hours}h {minutes}min"
            else:
                return f"{minutes}min"
        except:
            return "unknown"


def main():
    """Main entry point for CLI usage"""
    monitor = TrackingMonitor()

    if len(sys.argv) > 1:
        command = sys.argv[1]

        if command == "status":
            monitor.print_status()

        elif command == "check":
            action, reason = monitor.check_cycle()
            if action == "archive":
                print(f"⚠️  Archive needed: {reason}")
                sys.exit(1)
            else:
                print("✓ Continue working (no archive needed)")
                sys.exit(0)

        elif command == "violations":
            violations = monitor.enforce_rules()
            if violations:
                print(f"⚠️  {len(violations)} violation(s) detected:")
                for v in violations:
                    print(f"\n  {v['type'].upper()}")
                    print(f"  Files: {v.get('files', v.get('directories', []))}")
                    print(f"  Action: {v['action']}")
                sys.exit(1)
            else:
                print("✓ No violations detected")
                sys.exit(0)

        elif command == "archive":
            summary = sys.argv[2] if len(sys.argv) > 2 else None
            success = monitor.archive_cycle(summary)
            sys.exit(0 if success else 1)

        else:
            print(f"Unknown command: {command}")
            print("Usage: monitor.py [status|check|violations|archive]")
            sys.exit(1)

    else:
        # Default: print status
        monitor.print_status()


if __name__ == "__main__":
    main()
