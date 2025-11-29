#!/bin/bash
# Generate personalized student report for presentation
# Usage: ./student_report.sh <STUDENT_CODE>
#        ./student_report.sh --all
# Example: ./student_report.sh HV

set -e

CONTAINER="cosmelab/dna-barcoding-analysis:latest"
ALL_STUDENTS="BR HV JA JM JR KG MA TW WA WL"

if [ -z "$1" ]; then
    echo "Usage: ./student_report.sh <STUDENT_CODE>"
    echo "       ./student_report.sh --all     (generate all reports)"
    echo ""
    echo "Available student codes:"
    echo "  BR, HV, JA, JM, JR, KG, MA, TW, WA, WL"
    echo ""
    echo "Example: ./student_report.sh HV"
    exit 1
fi

# Create output directory
mkdir -p results/student_reports

# Start logging
LOG_FILE="results/student_reports/report_generation.log"
exec > >(tee "$LOG_FILE") 2>&1
echo "=== Student Report Generation ==="
echo "🖥️  Environment: Local (Docker)"
echo "📝 Log file: $LOG_FILE"
echo ""

if [ "$1" == "--all" ]; then
    # Show Rich banner
    docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
      $CONTAINER \
      python3 -c "
from rich.console import Console
from rich.panel import Panel
console = Console(force_terminal=True)
console.print(Panel('[bold cyan]GENERATING ALL STUDENT REPORTS[/]', border_style='cyan'))
"
    echo ""

    for STUDENT in $ALL_STUDENTS; do
        echo "📄 Generating report for: $STUDENT"
        docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace/modules/06_lab_data_analysis \
          $CONTAINER \
          python3 generate_student_report.py "$STUDENT"
        echo ""
    done

    # Show completion banner
    docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
      $CONTAINER \
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
    docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
      $CONTAINER \
      python3 -c "
from rich.console import Console
from rich.panel import Panel
console = Console(force_terminal=True)
console.print(Panel('[bold cyan]GENERATING STUDENT REPORT: $STUDENT_CODE[/]', border_style='cyan'))
"
    echo ""

    docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace/modules/06_lab_data_analysis \
      $CONTAINER \
      python3 generate_student_report.py "$STUDENT_CODE"

    echo ""
    # Show completion banner
    docker run --rm --entrypoint="" -v $(pwd):/workspace -w /workspace \
      $CONTAINER \
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
