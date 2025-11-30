#!/usr/bin/env python3
"""
Show Welcome Screen for DNA Barcoding Pipeline

This script displays a beautiful welcome message with:
- DNA Barcoding banner
- Workflow diagram
- Environment options

Usage:
    python3 modules/show_welcome.py [--step N] [--mode tutorial|analysis]
"""

import argparse
import sys

# Import our UI module
try:
    from pipeline_ui import (
        print_banner, print_workflow, print_environment_choice,
        print_step_banner, print_step_header, print_completion,
        print_team_battle, print_success, print_info
    )
except ImportError:
    # If run from different directory, try adding modules to path
    sys.path.insert(0, '/workspace/modules')
    sys.path.insert(0, 'modules')
    from pipeline_ui import (
        print_banner, print_workflow, print_environment_choice,
        print_step_banner, print_step_header, print_completion,
        print_team_battle, print_success, print_info
    )


def show_welcome(mode: str = "analysis"):
    """Show the welcome screen."""
    print_banner()
    print_workflow()

    if mode == "tutorial":
        print_info("Running TUTORIAL with test data...")
        print_info("This teaches you the workflow before analyzing your own data.\n")
    else:
        print_info("Running ANALYSIS on class data...")
        print_info("Make sure you've completed the tutorial first!\n")


def show_step(step_num: int, total_steps: int = 6):
    """Show a specific step banner."""
    step_titles = {
        1: "Quality Control",
        2: "Consensus Sequences",
        3: "Combine with References",
        4: "Alignment & Phylogenetic Tree",
        5: "Species Identification (BLAST)",
        6: "Lab Data Analysis"
    }
    title = step_titles.get(step_num, f"Step {step_num}")
    print_step_header(step_num, total_steps, title)
    print_step_banner(step_num)


def show_completion():
    """Show completion message."""
    print_completion()


def show_teams():
    """Show team battle banner."""
    print_team_battle()


def main():
    parser = argparse.ArgumentParser(description="Show pipeline UI elements")
    parser.add_argument("--welcome", action="store_true", help="Show welcome screen")
    parser.add_argument("--step", type=int, help="Show step N banner (1-6)")
    parser.add_argument("--complete", action="store_true", help="Show completion banner")
    parser.add_argument("--teams", action="store_true", help="Show team battle")
    parser.add_argument("--mode", choices=["tutorial", "analysis"], default="analysis",
                        help="Mode for welcome message")
    parser.add_argument("--workflow", action="store_true", help="Show workflow diagram only")

    args = parser.parse_args()

    if args.welcome:
        show_welcome(args.mode)
    elif args.step:
        show_step(args.step)
    elif args.complete:
        show_completion()
    elif args.teams:
        show_teams()
    elif args.workflow:
        print_workflow()
    else:
        # Default: show everything
        show_welcome(args.mode)


if __name__ == "__main__":
    main()
