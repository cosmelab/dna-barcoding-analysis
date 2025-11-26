#!/usr/bin/env python3
"""
Interactive Assignment Question Script

Asks students questions about their DNA barcoding analysis results.
Saves answers to answers.json for automated grading.

Run this AFTER completing ./tutorial.sh and ./run-analysis.sh
"""

import json
import sys
from pathlib import Path
from typing import Dict, List


class AssignmentQuestions:
    """Interactive assignment questions."""

    def __init__(self):
        self.answers = {}
        self.results_dir = Path("results/my_analysis")

    def print_header(self):
        """Print welcome message."""
        print("\n" + "=" * 70)
        print("DNA BARCODING ASSIGNMENT - INTERACTIVE QUESTIONS")
        print("=" * 70)
        print("\nThis script will ask you questions about your analysis results.")
        print("Your answers will be saved to answers.json for grading.\n")
        print("Make sure you've completed:")
        print("  1. ./tutorial.sh")
        print("  2. ./run-analysis.sh\n")

        input("Press ENTER to start...")
        print()

    def check_analysis_complete(self) -> bool:
        """Check if analysis has been run."""
        required_files = [
            self.results_dir / "qc" / "qc_report.html",
            self.results_dir / "consensus" / "consensus_sequences.fasta",
            self.results_dir / "phylogeny" / "tree.png",
            self.results_dir / "blast" / "identification_report.html"
        ]

        missing = [f for f in required_files if not f.exists()]

        if missing:
            print("❌ ERROR: Analysis not complete. Missing files:")
            for f in missing:
                print(f"   - {f}")
            print("\nPlease run ./run-analysis.sh first!")
            return False

        return True

    def ask_multiple_choice(self, question: str, options: List[str]) -> str:
        """Ask multiple choice question."""
        print(f"\n{question}\n")
        for i, option in enumerate(options, 1):
            print(f"  {i}. {option}")

        while True:
            try:
                choice = input("\nYour answer (enter number): ").strip()
                idx = int(choice) - 1
                if 0 <= idx < len(options):
                    return options[idx]
                else:
                    print(f"Please enter a number between 1 and {len(options)}")
            except ValueError:
                print("Please enter a valid number")

    def ask_open_ended(self, question: str, hint: str = "") -> str:
        """Ask open-ended question."""
        print(f"\n{question}")
        if hint:
            print(f"   Hint: {hint}")
        print()

        answer = input("Your answer: ").strip()
        while not answer:
            print("Please provide an answer.")
            answer = input("Your answer: ").strip()

        return answer

    def ask_yes_no(self, question: str) -> bool:
        """Ask yes/no question."""
        print(f"\n{question}")

        while True:
            answer = input("Your answer (yes/no): ").strip().lower()
            if answer in ['yes', 'y']:
                return True
            elif answer in ['no', 'n']:
                return False
            else:
                print("Please answer 'yes' or 'no'")

    def section_results_table(self):
        """Part 2: Results Table questions."""
        print("\n" + "=" * 70)
        print("PART 2: SPECIES IDENTIFICATION RESULTS")
        print("=" * 70)
        print("\nOpen the BLAST report: results/my_analysis/blast/identification_report.html")
        print("Look at the table showing species identifications.\n")

        input("Press ENTER when you've opened the report...")

        # Ask how many samples have results
        num_samples = self.ask_open_ended(
            "How many samples have BLAST results (consensus sequences)?",
            hint="Count the rows in the table"
        )
        self.answers['num_samples'] = num_samples

        # Collect species data for each sample
        self.answers['species_table'] = []

        try:
            n = int(num_samples)
            for i in range(n):
                print(f"\n--- Sample {i+1} of {n} ---")

                sample_name = self.ask_open_ended(f"Sample name (e.g., AT-HV1):")

                species = self.ask_open_ended(
                    f"Top BLAST hit species for {sample_name}:",
                    hint="Format: Genus species (e.g., Culex quinquefasciatus)"
                )

                percent_id = self.ask_open_ended(
                    f"% Identity for {sample_name}:",
                    hint="Format: 99.5 (just the number)"
                )

                common_name = self.ask_open_ended(
                    f"Common name for {sample_name} (optional, or type 'unknown'):"
                )

                self.answers['species_table'].append({
                    'sample': sample_name,
                    'species': species,
                    'percent_identity': percent_id,
                    'common_name': common_name if common_name.lower() != 'unknown' else ''
                })

        except ValueError:
            print("Invalid number of samples. Skipping table entry.")

    def section_quality_control(self):
        """Part 3: Question 1 - Quality Control."""
        print("\n" + "=" * 70)
        print("PART 3: QUESTION 1 - QUALITY CONTROL")
        print("=" * 70)
        print("\nOpen the QC report: results/my_analysis/qc/qc_report.html\n")

        input("Press ENTER when you've opened the report...")

        # Q1a: Sequences passed
        q1a = self.ask_open_ended(
            "Q1(a): How many of the 30 class sequences (.ab1 files) passed QC?",
            hint="Look at the summary at the top of the report"
        )
        self.answers['q1a_sequences_passed'] = q1a

        # Q1b: Pairs passed
        q1b = self.ask_open_ended(
            "Q1(b): How many samples had BOTH forward AND reverse reads pass?",
            hint="Count samples where both F and R are in the 'Passed' section"
        )
        self.answers['q1b_pairs_passed'] = q1b

        # Q1c: Why F+R important
        q1c = self.ask_open_ended(
            "Q1(c): Why is it important to have both F and R reads?",
            hint="Think about coverage, accuracy, error detection"
        )
        self.answers['q1c_explanation'] = q1c

    def section_phylogenetic_tree(self):
        """Part 3: Question 2 - Phylogenetic Tree."""
        print("\n" + "=" * 70)
        print("PART 3: QUESTION 2 - PHYLOGENETIC TREE")
        print("=" * 70)
        print("\nOpen the tree: results/my_analysis/phylogeny/tree.png\n")

        input("Press ENTER when you've opened the tree...")

        # Q2a: Clustering pattern
        q2a_options = [
            "Class samples cluster together in one group",
            "Class samples are spread across different parts of the tree",
            "Class samples form multiple separate clusters"
        ]
        q2a = self.ask_multiple_choice(
            "Q2(a): Do the class samples cluster together, or are they spread across different parts of the tree?",
            q2a_options
        )
        self.answers['q2a_clustering'] = q2a

        # Q2b: Related species
        q2b = self.ask_open_ended(
            "Q2(b): Which reference species are the class samples most closely related to?",
            hint="Look at which reference species are nearest to your samples on the tree"
        )
        self.answers['q2b_related_species'] = q2b

        # Q2c: Diversity interpretation
        q2c = self.ask_open_ended(
            "Q2(c): What does this tell you about mosquito diversity in the class sampling locations?",
            hint="Multiple species = high diversity, one species = low diversity"
        )
        self.answers['q2c_diversity'] = q2c

    def section_species_identification(self):
        """Part 3: Question 3 - Species Identification."""
        print("\n" + "=" * 70)
        print("PART 3: QUESTION 3 - SPECIES IDENTIFICATION")
        print("=" * 70)

        # Q3a: List species
        q3a = self.ask_open_ended(
            "Q3(a): What mosquito species did the class identify? List all unique species found.",
            hint="Format: Genus species, Genus species (comma separated)"
        )
        self.answers['q3a_species_list'] = q3a

        # Q3b: BLAST vs tree agreement
        q3b = self.ask_yes_no(
            "Q3(b): Do the BLAST results (% identity) agree with where samples clustered on the tree?"
        )
        self.answers['q3b_blast_tree_agree'] = "yes" if q3b else "no"

        # Q3c: Southern California
        q3c = self.ask_yes_no(
            "Q3(c): Are these species known to occur in Southern California?"
        )
        self.answers['q3c_socal_species'] = "yes" if q3c else "no"

        if not q3c:
            q3c_explain = self.ask_open_ended(
                "Please explain which species are NOT found in Southern California:"
            )
            self.answers['q3c_explanation'] = q3c_explain

        # Q3d: Confidence
        q3d_options = [
            "Very confident (>99% identity, matches known species)",
            "Confident (>98% identity, matches expected species)",
            "Somewhat confident (95-98% identity)",
            "Not confident (<95% identity or unexpected species)"
        ]
        q3d = self.ask_multiple_choice(
            "Q3(d): How confident are you in the species identifications?",
            q3d_options
        )
        self.answers['q3d_confidence'] = q3d

    def save_answers(self):
        """Save answers to JSON file."""
        output_file = Path("answers.json")

        with open(output_file, 'w') as f:
            json.dump(self.answers, f, indent=2)

        print("\n" + "=" * 70)
        print("✓ ANSWERS SAVED!")
        print("=" * 70)
        print(f"\nYour answers have been saved to: {output_file}")
        print("\nNext steps:")
        print("  1. git add answers.json")
        print("  2. git commit -m 'Complete assignment answers'")
        print("  3. git push origin main")
        print("\nThe auto-grader will check your answers when you push.\n")

    def run(self):
        """Run the interactive questionnaire."""
        self.print_header()

        if not self.check_analysis_complete():
            sys.exit(1)

        # Run all question sections
        self.section_results_table()
        self.section_quality_control()
        self.section_phylogenetic_tree()
        self.section_species_identification()

        # Save results
        self.save_answers()

        print("Good luck! 🧬🦟\n")


def main():
    """Main entry point."""
    questionnaire = AssignmentQuestions()
    try:
        questionnaire.run()
    except KeyboardInterrupt:
        print("\n\n❌ Cancelled. Your answers were not saved.")
        sys.exit(1)
    except Exception as e:
        print(f"\n\n❌ Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
