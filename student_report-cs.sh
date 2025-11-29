#!/bin/bash
# Generate personalized student report for presentation (Codespaces Version)
# Usage: ./student_report-cs.sh <STUDENT_CODE>
#        ./student_report-cs.sh --all
# Example: ./student_report-cs.sh HV

set -e

ALL_STUDENTS="BR HV JA JM JR KG MA TW WA WL"

if [ -z "$1" ]; then
    echo "Usage: ./student_report-cs.sh <STUDENT_CODE>"
    echo "       ./student_report-cs.sh --all     (generate all reports)"
    echo ""
    echo "Available student codes:"
    echo "  BR, HV, JA, JM, JR, KG, MA, TW, WA, WL"
    echo ""
    echo "Example: ./student_report-cs.sh HV"
    exit 1
fi

# Create output directory
mkdir -p results/student_reports

# Start logging
LOG_FILE="results/student_reports/report_generation.log"
exec > >(tee "$LOG_FILE") 2>&1
echo "=== Student Report Generation ==="
echo "☁️  Environment: GitHub Codespaces"
echo "📝 Log file: $LOG_FILE"
echo ""

if [ "$1" == "--all" ]; then
    # Show Rich banner
    python3 -c "
from rich.console import Console
from rich.panel import Panel
console = Console(force_terminal=True)
console.print(Panel('[bold cyan]GENERATING ALL STUDENT REPORTS[/]', border_style='cyan'))
"
    echo ""

    for STUDENT in $ALL_STUDENTS; do
        echo "📄 Generating report for: $STUDENT"
        cd modules/06_lab_data_analysis
        python3 generate_student_report.py "$STUDENT"
        cd ../..
        echo ""
    done

    # Show completion banner
    python3 -c "
from rich.console import Console
from rich.panel import Panel
console = Console(force_terminal=True)
console.print(Panel('[bold green]✓ ALL REPORTS GENERATED![/]', border_style='green'))
"
    echo ""
    echo "Reports saved to: results/student_reports/"
    echo ""
    ls -la results/student_reports/*.html
    echo ""
else
    STUDENT_CODE="$1"

    # Show Rich banner
    python3 -c "
from rich.console import Console
from rich.panel import Panel
console = Console(force_terminal=True)
console.print(Panel('[bold cyan]GENERATING STUDENT REPORT: $STUDENT_CODE[/]', border_style='cyan'))
"
    echo ""

    cd modules/06_lab_data_analysis
    python3 generate_student_report.py "$STUDENT_CODE"
    cd ../..

    echo ""
    # Show completion banner
    python3 -c "
from rich.console import Console
from rich.panel import Panel
console = Console(force_terminal=True)
console.print(Panel('[bold green]✓ REPORT GENERATED![/]', border_style='green'))
"
    echo ""
    echo "📂 Open: results/student_reports/${STUDENT_CODE}_report.html"
    echo ""
fi

echo "=== Report Generation Completed ==="
