#!/usr/bin/env python3
"""
Interactive Assignment - DNA Barcoding Analysis

A LEARNING-focused assignment that validates your understanding.
Uses the class dataset - everyone has the same correct answers!

Features:
- Retry mechanism (max 3 attempts per question)
- Progressive hints that get more specific
- Immediate feedback with explanations
- Beautiful Rich terminal UI
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    from rich.console import Console
    from rich.panel import Panel
    from rich.table import Table
    from rich.text import Text
    from rich.prompt import Prompt, IntPrompt
    from rich import box
    RICH_AVAILABLE = True
except ImportError:
    RICH_AVAILABLE = False

# Initialize console
console = Console() if RICH_AVAILABLE else None


# ============================================================================
# ANSWER KEY - Everyone analyzes the same class data
# ============================================================================
ANSWER_KEY = {
    'q1a_sequences_passed': {
        'correct': ['12'],
        'hints': [
            "Look at the QC report summary at the top",
            "Count the sequences in the 'Passed QC' section",
            "The answer is between 10-15"
        ],
        'explanation': "12 out of 30 sequences passed our quality thresholds (Q20 average, 400+ good bases)."
    },
    'q1b_pairs_passed': {
        'correct': ['4'],
        'hints': [
            "A pair means BOTH forward (F) AND reverse (R) passed",
            "Look for samples where you see both F and R in the passed list",
            "Count: AT-HV1, AT-HV3, AT-JM2, AT-WL2"
        ],
        'explanation': "Only 4 samples had complete F+R pairs. This is why consensus sequences are so valuable - they require both reads!"
    },
    'num_samples': {
        'correct': ['4'],
        'hints': [
            "This should match the number of pairs that passed QC",
            "Count the rows in the BLAST results table",
            "Each consensus sequence = one BLAST result"
        ],
        'explanation': "4 consensus sequences were created from the 4 complete F+R pairs."
    }
}


def print_header():
    """Print beautiful welcome header."""
    if RICH_AVAILABLE:
        console.print()
        console.print(Panel.fit(
            "[bold cyan]DNA BARCODING ASSIGNMENT[/]\n"
            "[dim]Interactive Learning Experience[/]",
            border_style="cyan",
            box=box.DOUBLE
        ))
        console.print()
        console.print("[green]This is a LEARNING tool![/] You can retry questions if you get them wrong.")
        console.print("Everyone analyzes the [bold]same class data[/], so we can validate your answers.\n")
        console.print("[dim]Tip: Keep your HTML reports open in a browser tab.[/]\n")
    else:
        print("\n" + "=" * 60)
        print("DNA BARCODING ASSIGNMENT")
        print("Interactive Learning Experience")
        print("=" * 60)
        print("\nThis is a LEARNING tool! You can retry questions.")
        print()


def print_section(title: str, subtitle: str = ""):
    """Print section header."""
    if RICH_AVAILABLE:
        console.print()
        console.rule(f"[bold yellow]{title}[/]")
        if subtitle:
            console.print(f"[dim]{subtitle}[/]")
        console.print()
    else:
        print("\n" + "=" * 60)
        print(title)
        if subtitle:
            print(subtitle)
        print("=" * 60 + "\n")


def print_success(message: str):
    """Print success message."""
    if RICH_AVAILABLE:
        console.print(f"[bold green]Correct![/] {message}")
    else:
        print(f"Correct! {message}")


def print_error(message: str, hint: str = ""):
    """Print error message with hint."""
    if RICH_AVAILABLE:
        console.print(f"[bold red]Not quite.[/] {message}")
        if hint:
            console.print(f"[yellow]Hint:[/] {hint}")
    else:
        print(f"Not quite. {message}")
        if hint:
            print(f"Hint: {hint}")


def print_final_answer(correct_answer: str, explanation: str):
    """Print the correct answer after max attempts."""
    if RICH_AVAILABLE:
        console.print(Panel(
            f"[bold]Correct answer:[/] {correct_answer}\n\n"
            f"[dim]{explanation}[/]",
            title="[yellow]Learning moment[/]",
            border_style="yellow"
        ))
    else:
        print(f"\nCorrect answer: {correct_answer}")
        print(f"Explanation: {explanation}\n")


def ask_validated_question(
    question: str,
    key: str,
    max_attempts: int = 3
) -> Tuple[str, int, bool]:
    """
    Ask a question with validation and retry.

    Returns: (answer, attempts, was_correct)
    """
    config = ANSWER_KEY.get(key, {})
    correct_answers = config.get('correct', [])
    hints = config.get('hints', ["Check the report again"])
    explanation = config.get('explanation', "")

    for attempt in range(max_attempts):
        # Show question
        if RICH_AVAILABLE:
            if attempt > 0:
                console.print(f"[dim]Attempt {attempt + 1} of {max_attempts}[/]")
            console.print(f"[bold]{question}[/]")
            answer = Prompt.ask("[cyan]Your answer[/]")
        else:
            if attempt > 0:
                print(f"(Attempt {attempt + 1} of {max_attempts})")
            print(question)
            answer = input("Your answer: ").strip()

        # Check if correct
        if answer.strip() in correct_answers:
            print_success(explanation)
            return answer, attempt + 1, True

        # Wrong answer
        if attempt < max_attempts - 1:
            hint = hints[min(attempt, len(hints) - 1)]
            print_error("Try again!", hint)
            print()
        else:
            # Max attempts reached - show answer
            print_final_answer(correct_answers[0], explanation)
            return answer, max_attempts, False

    return answer, max_attempts, False


def ask_open_question(question: str, hint: str = "") -> str:
    """Ask open-ended question (no validation)."""
    if RICH_AVAILABLE:
        console.print(f"[bold]{question}[/]")
        if hint:
            console.print(f"[dim]Hint: {hint}[/]")
        return Prompt.ask("[cyan]Your answer[/]")
    else:
        print(question)
        if hint:
            print(f"Hint: {hint}")
        return input("Your answer: ").strip()


def ask_multiple_choice(question: str, options: List[str]) -> str:
    """Ask multiple choice question."""
    if RICH_AVAILABLE:
        console.print(f"[bold]{question}[/]\n")
        for i, opt in enumerate(options, 1):
            console.print(f"  [cyan]{i}.[/] {opt}")

        while True:
            try:
                choice = IntPrompt.ask("\n[cyan]Enter number[/]")
                if 1 <= choice <= len(options):
                    return options[choice - 1]
                console.print(f"[red]Please enter 1-{len(options)}[/]")
            except:
                console.print("[red]Please enter a number[/]")
    else:
        print(f"{question}\n")
        for i, opt in enumerate(options, 1):
            print(f"  {i}. {opt}")

        while True:
            try:
                choice = int(input("\nEnter number: "))
                if 1 <= choice <= len(options):
                    return options[choice - 1]
            except ValueError:
                pass
            print(f"Please enter 1-{len(options)}")


def check_analysis_complete() -> bool:
    """Check if analysis has been run."""
    results_dir = Path("results/my_analysis")
    required = [
        results_dir / "01_qc" / "qc_report.html",
        results_dir / "05_blast" / "identification_report.html",
        results_dir / "04_phylogeny" / "tree.png"
    ]

    missing = [f for f in required if not f.exists()]

    if missing:
        if RICH_AVAILABLE:
            console.print("[bold red]Analysis not complete![/]")
            console.print("\nMissing files:")
            for f in missing:
                console.print(f"  [red]x[/] {f}")
            console.print("\n[yellow]Run ./run-analysis.sh first![/]")
        else:
            print("Analysis not complete!")
            for f in missing:
                print(f"  Missing: {f}")
            print("\nRun ./run-analysis.sh first!")
        return False

    return True


def run_assignment():
    """Run the interactive assignment."""
    answers = {
        'metadata': {
            'version': '2.0',
            'learning_mode': True
        }
    }

    print_header()

    if not check_analysis_complete():
        sys.exit(1)

    if RICH_AVAILABLE:
        Prompt.ask("[dim]Press Enter to begin[/]")
    else:
        input("Press Enter to begin...")

    # ========================================================================
    # SECTION 1: Quality Control
    # ========================================================================
    print_section(
        "PART 1: Quality Control",
        "Open: results/my_analysis/01_qc/qc_report.html"
    )

    # Q1a: Sequences passed (validated)
    answer, attempts, correct = ask_validated_question(
        "How many sequences PASSED quality control?",
        'q1a_sequences_passed'
    )
    answers['q1a_sequences_passed'] = answer
    answers['q1a_attempts'] = attempts
    answers['q1a_correct'] = correct

    print()

    # Q1b: Pairs passed (validated)
    answer, attempts, correct = ask_validated_question(
        "How many samples had BOTH forward AND reverse reads pass?",
        'q1b_pairs_passed'
    )
    answers['q1b_pairs_passed'] = answer
    answers['q1b_attempts'] = attempts
    answers['q1b_correct'] = correct

    print()

    # Q1c: Why F+R important (open)
    answers['q1c_explanation'] = ask_open_question(
        "Why is it important to have both forward and reverse reads?",
        "Think about coverage, accuracy, and error correction"
    )

    # ========================================================================
    # SECTION 2: Species Identification
    # ========================================================================
    print_section(
        "PART 2: Species Identification",
        "Open: results/my_analysis/05_blast/identification_report.html"
    )

    # Number of samples (validated)
    answer, attempts, correct = ask_validated_question(
        "How many consensus sequences have BLAST results?",
        'num_samples'
    )
    answers['num_samples'] = answer
    answers['num_samples_attempts'] = attempts

    # Collect species table
    answers['species_table'] = []

    try:
        n = int(answer)
        if RICH_AVAILABLE:
            console.print(f"\n[dim]Enter details for each of the {n} samples:[/]\n")
        else:
            print(f"\nEnter details for each of the {n} samples:\n")

        for i in range(n):
            if RICH_AVAILABLE:
                console.rule(f"[dim]Sample {i+1} of {n}[/]")
            else:
                print(f"--- Sample {i+1} of {n} ---")

            sample = ask_open_question("Sample name (e.g., AT-HV1):")
            species = ask_open_question(
                "Top BLAST hit species:",
                "Format: Genus species (e.g., Aedes albopictus)"
            )
            percent = ask_open_question(
                "Percent identity:",
                "Just the number (e.g., 99.5)"
            )

            answers['species_table'].append({
                'sample': sample,
                'species': species,
                'percent_identity': percent
            })
            print()
    except ValueError:
        pass

    # ========================================================================
    # SECTION 3: Phylogenetic Analysis
    # ========================================================================
    print_section(
        "PART 3: Phylogenetic Analysis",
        "Open: results/my_analysis/04_phylogeny/tree.png"
    )

    # Q2a: Clustering pattern
    answers['q2a_clustering'] = ask_multiple_choice(
        "How are the class samples distributed on the tree?",
        [
            "Clustered together in one group",
            "Spread across different parts (multiple species)",
            "Cannot determine from the tree"
        ]
    )

    print()

    # Q2b: Related species
    answers['q2b_related_species'] = ask_open_question(
        "Which reference species are the class samples most closely related to?",
        "Look at which references cluster near your samples"
    )

    print()

    # Q2c: Diversity
    answers['q2c_diversity'] = ask_open_question(
        "What does this tell you about mosquito diversity in the sampling locations?",
        "Multiple species = higher diversity"
    )

    # ========================================================================
    # SECTION 4: Interpretation
    # ========================================================================
    print_section("PART 4: Final Questions")

    # Q3a: All species
    answers['q3a_species_list'] = ask_open_question(
        "List ALL unique species identified by the class:",
        "Comma-separated: Species1, Species2"
    )

    print()

    # Q3b: Confidence
    answers['q3b_confidence'] = ask_multiple_choice(
        "How confident are you in these identifications?",
        [
            "Very confident (>99% identity)",
            "Confident (97-99% identity)",
            "Somewhat confident (95-97% identity)",
            "Not confident (<95% identity)"
        ]
    )

    # ========================================================================
    # SAVE ANSWERS
    # ========================================================================
    output_dir = Path("submission")
    output_dir.mkdir(exist_ok=True)
    output_file = output_dir / "answers.json"
    with open(output_file, 'w') as f:
        json.dump(answers, f, indent=2)

    if RICH_AVAILABLE:
        console.print()
        console.print(Panel.fit(
            "[bold green]Assignment Complete![/]\n\n"
            f"Answers saved to: [cyan]{output_file}[/]\n\n"
            "[dim]Next steps:[/]\n"
            "  1. git add submission/answers.json\n"
            "  2. git commit -m 'Complete assignment'\n"
            "  3. git push",
            border_style="green",
            box=box.DOUBLE
        ))
    else:
        print("\n" + "=" * 60)
        print("Assignment Complete!")
        print("=" * 60)
        print(f"\nAnswers saved to: {output_file}")
        print("\nNext steps:")
        print("  1. git add submission/answers.json")
        print("  2. git commit -m 'Complete assignment'")
        print("  3. git push")

    print()


def main():
    """Main entry point."""
    try:
        run_assignment()
    except KeyboardInterrupt:
        if RICH_AVAILABLE:
            console.print("\n[yellow]Cancelled. Answers not saved.[/]")
        else:
            print("\nCancelled. Answers not saved.")
        sys.exit(1)


if __name__ == "__main__":
    main()
