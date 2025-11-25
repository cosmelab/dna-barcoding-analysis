#!/usr/bin/env python3
"""
Shared utility functions for DNA barcoding analysis modules
Provides consistent terminal output formatting with rich library and browser opening
"""

from pathlib import Path
import webbrowser


def print_header(text):
    """Print a formatted header"""
    width = 70
    print("\n" + "=" * width)
    print(f"  {text}")
    print("=" * width)


def print_step(step_num, total_steps, description):
    """Print a step indicator"""
    print(f"\nStep {step_num}/{total_steps}: {description}")
    print("-" * 70)


def print_success(message):
    """Print a success message"""
    print(f"✓ {message}")


def print_info(message, indent=True):
    """Print an info message"""
    prefix = "  " if indent else ""
    print(f"{prefix}{message}")


def print_error(message):
    """Print an error message"""
    print(f"✗ ERROR: {message}")


def print_warning(message):
    """Print a warning message"""
    print(f"⚠ WARNING: {message}")


def open_in_browser(file_path):
    """
    Open HTML file in default web browser (cross-platform)

    Works on Mac, Windows, and Linux. Handles both absolute
    and relative paths, and converts to proper file:// URLs.

    Args:
        file_path: Path to HTML file (str or Path object)

    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # Convert to absolute path and use file:// URL
        abs_path = Path(file_path).resolve()

        # On Windows, webbrowser needs forward slashes even though Path uses backslashes
        url = f'file://{abs_path.as_posix()}'

        webbrowser.open(url)
        return True
    except Exception as e:
        print(f"Could not auto-open browser: {e}")
        return False


def print_next_steps(steps):
    """
    Print a formatted "NEXT STEPS" section

    Args:
        steps: List of step dictionaries with 'title' and 'details'
               Example: [{'title': 'View results', 'details': '/path/to/file.html'}]
    """
    print("\n" + "=" * 70)
    print("  NEXT STEPS:")
    print("=" * 70)

    for i, step in enumerate(steps, 1):
        print_info(f"{i}. {step['title']}", indent=False)
        if isinstance(step['details'], list):
            for detail in step['details']:
                print_info(f"   {detail}", indent=False)
        else:
            print_info(f"   {step['details']}", indent=False)
        if i < len(steps):
            print_info("", indent=False)

    print("=" * 70 + "\n")
