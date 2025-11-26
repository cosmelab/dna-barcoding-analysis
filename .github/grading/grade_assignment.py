#!/usr/bin/env python3
"""
Automated Grading Script for DNA Barcoding Assignment

Parses student assignment.md and grades based on:
- Results table (species IDs, % identity)
- QC statistics
- Written answers (keyword matching)

Since all students analyze the same class dataset, answers should be identical.
"""

import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple


class AssignmentGrader:
    """Grade DNA barcoding assignment answers."""

    def __init__(self, assignment_path: str, answer_key_path: str):
        self.assignment_path = Path(assignment_path)
        self.answer_key_path = Path(answer_key_path)
        self.total_points = 0
        self.earned_points = 0
        self.feedback = []

    def load_assignment(self) -> str:
        """Load assignment markdown file."""
        with open(self.assignment_path) as f:
            return f.read()

    def load_answer_key(self) -> Dict:
        """Load answer key with expected responses."""
        with open(self.answer_key_path) as f:
            import json
            return json.load(f)

    def extract_results_table(self, content: str) -> List[Dict]:
        """Extract species identification table from markdown."""
        # Find the table in Part 2
        table_pattern = r'\| Sample \| Species Identified \| % Identity \| Common Name \|(.*?)---'
        match = re.search(table_pattern, content, re.DOTALL)

        if not match:
            return []

        table_content = match.group(1)
        rows = []

        # Parse each row
        for line in table_content.strip().split('\n'):
            if line.strip().startswith('|'):
                cells = [cell.strip() for cell in line.split('|')[1:-1]]
                if len(cells) == 4 and cells[0]:  # Skip header and empty rows
                    rows.append({
                        'sample': cells[0],
                        'species': cells[1],
                        'percent_identity': cells[2],
                        'common_name': cells[3]
                    })

        return rows

    def extract_question_answers(self, content: str) -> Dict[str, str]:
        """Extract answers from Question 1, 2, 3."""
        answers = {}

        # Question 1
        q1_pattern = r'### Question 1:.*?\*\*Your answer:\*\*\s*```(.*?)```'
        q1_match = re.search(q1_pattern, content, re.DOTALL)
        if q1_match:
            answers['question1'] = q1_match.group(1).strip()

        # Question 2
        q2_pattern = r'### Question 2:.*?\*\*Your answer:\*\*\s*```(.*?)```'
        q2_match = re.search(q2_pattern, content, re.DOTALL)
        if q2_match:
            answers['question2'] = q2_match.group(1).strip()

        # Question 3
        q3_pattern = r'### Question 3:.*?\*\*Your answer:\*\*\s*```(.*?)```'
        q3_match = re.search(q3_pattern, content, re.DOTALL)
        if q3_match:
            answers['question3'] = q3_match.group(1).strip()

        return answers

    def grade_results_table(self, table: List[Dict], answer_key: Dict) -> int:
        """Grade the species identification table (20 points)."""
        points = 0
        max_points = 20

        expected = answer_key['results_table']

        if len(table) == 0:
            self.feedback.append("❌ Results table is empty (0/20 points)")
            return 0

        # Check each species
        for i, expected_row in enumerate(expected):
            if i >= len(table):
                self.feedback.append(f"❌ Missing row for {expected_row['sample']}")
                continue

            student_row = table[i]
            row_correct = True

            # Check species name (case-insensitive, handle italics)
            student_species = student_row['species'].replace('*', '').strip().lower()
            expected_species = expected_row['species'].lower()

            if student_species != expected_species:
                self.feedback.append(
                    f"❌ {expected_row['sample']}: Wrong species. "
                    f"Expected {expected_row['species']}, got {student_row['species']}"
                )
                row_correct = False

            # Check % identity (allow ±1% tolerance)
            try:
                student_pct = float(student_row['percent_identity'].replace('%', ''))
                expected_pct = float(expected_row['percent_identity'])

                if abs(student_pct - expected_pct) > 1.0:
                    self.feedback.append(
                        f"❌ {expected_row['sample']}: Wrong % identity. "
                        f"Expected ~{expected_pct}%, got {student_pct}%"
                    )
                    row_correct = False
            except ValueError:
                self.feedback.append(f"❌ {expected_row['sample']}: Invalid % identity format")
                row_correct = False

            if row_correct:
                points += max_points / len(expected)

        self.feedback.append(f"✓ Results table: {points:.1f}/{max_points} points")
        return int(points)

    def grade_question1(self, answer: str, answer_key: Dict) -> int:
        """Grade Question 1 - Quality Control (5 points)."""
        points = 0
        max_points = 5

        expected = answer_key['question1']

        # Part a: Number of sequences passed (2 points)
        if str(expected['sequences_passed']) in answer:
            points += 2
            self.feedback.append("✓ Q1(a): Correct number of sequences passed QC")
        else:
            self.feedback.append(
                f"❌ Q1(a): Wrong answer. Expected {expected['sequences_passed']}"
            )

        # Part b: Pairs with both F and R (2 points)
        if str(expected['pairs_passed']) in answer:
            points += 2
            self.feedback.append("✓ Q1(b): Correct number of pairs")
        else:
            self.feedback.append(
                f"❌ Q1(b): Wrong answer. Expected {expected['pairs_passed']}"
            )

        # Part c: Conceptual understanding (1 point)
        keywords = ['coverage', 'accuracy', 'error', 'consensus', 'quality']
        if any(kw in answer.lower() for kw in keywords):
            points += 1
            self.feedback.append("✓ Q1(c): Good explanation of F+R importance")
        else:
            self.feedback.append("❌ Q1(c): Missing key concepts (coverage, accuracy, consensus)")

        self.feedback.append(f"✓ Question 1: {points}/{max_points} points")
        return points

    def grade_question2(self, answer: str, answer_key: Dict) -> int:
        """Grade Question 2 - Phylogenetic Tree (7 points)."""
        points = 0
        max_points = 7

        expected = answer_key['question2']

        # Part a: Clustering pattern (2 points)
        if any(kw in answer.lower() for kw in expected['clustering_keywords']):
            points += 2
            self.feedback.append("✓ Q2(a): Correctly described clustering pattern")
        else:
            self.feedback.append("❌ Q2(a): Missing description of sample clustering")

        # Part b: Reference species (3 points)
        mentioned_species = sum(1 for sp in expected['related_species'] if sp.lower() in answer.lower())
        if mentioned_species >= len(expected['related_species']) * 0.7:  # 70% of species
            points += 3
            self.feedback.append("✓ Q2(b): Correctly identified related species")
        else:
            self.feedback.append("❌ Q2(b): Missing or incorrect reference species")

        # Part c: Diversity interpretation (2 points)
        if any(kw in answer.lower() for kw in expected['diversity_keywords']):
            points += 2
            self.feedback.append("✓ Q2(c): Good interpretation of mosquito diversity")
        else:
            self.feedback.append("❌ Q2(c): Weak diversity interpretation")

        self.feedback.append(f"✓ Question 2: {points}/{max_points} points")
        return points

    def grade_question3(self, answer: str, answer_key: Dict) -> int:
        """Grade Question 3 - Species Identification (8 points)."""
        points = 0
        max_points = 8

        expected = answer_key['question3']

        # Part a: List species (2 points)
        mentioned_species = sum(1 for sp in expected['species_list'] if sp.lower() in answer.lower())
        if mentioned_species >= len(expected['species_list']) * 0.8:  # 80% of species
            points += 2
            self.feedback.append("✓ Q3(a): Correctly listed identified species")
        else:
            self.feedback.append("❌ Q3(a): Missing or incorrect species list")

        # Part b: BLAST vs tree agreement (2 points)
        if any(kw in answer.lower() for kw in ['agree', 'match', 'consistent', 'yes']):
            points += 2
            self.feedback.append("✓ Q3(b): Correctly noted BLAST-tree agreement")
        else:
            self.feedback.append("❌ Q3(b): Missing BLAST-tree comparison")

        # Part c: Southern California occurrence (2 points)
        if 'yes' in answer.lower() or 'southern california' in answer.lower():
            points += 2
            self.feedback.append("✓ Q3(c): Confirmed SoCal occurrence")
        else:
            self.feedback.append("❌ Q3(c): Missing geographic verification")

        # Part d: Confidence assessment (2 points)
        if any(kw in answer.lower() for kw in ['high', 'confident', '98', '99', '100']):
            points += 2
            self.feedback.append("✓ Q3(d): Good confidence assessment")
        else:
            self.feedback.append("❌ Q3(d): Missing confidence discussion")

        self.feedback.append(f"✓ Question 3: {points}/{max_points} points")
        return points

    def grade(self) -> Tuple[int, int, List[str]]:
        """Grade the assignment and return (earned, total, feedback)."""
        content = self.load_assignment()
        answer_key = self.load_answer_key()

        # Extract student answers
        table = self.extract_results_table(content)
        questions = self.extract_question_answers(content)

        # Grade each part
        self.feedback.append("\n## GRADING RESULTS\n")

        # Part 2: Results Table (20 points)
        self.earned_points += self.grade_results_table(table, answer_key)

        # Part 3: Questions
        if 'question1' in questions:
            self.earned_points += self.grade_question1(questions['question1'], answer_key)
        else:
            self.feedback.append("❌ Question 1: No answer found (0/5 points)")

        if 'question2' in questions:
            self.earned_points += self.grade_question2(questions['question2'], answer_key)
        else:
            self.feedback.append("❌ Question 2: No answer found (0/7 points)")

        if 'question3' in questions:
            self.earned_points += self.grade_question3(questions['question3'], answer_key)
        else:
            self.feedback.append("❌ Question 3: No answer found (0/8 points)")

        self.total_points = 40  # Part 2 (20) + Part 3 (20)

        # Final summary
        percentage = (self.earned_points / self.total_points) * 100
        self.feedback.append(f"\n## FINAL SCORE: {self.earned_points}/{self.total_points} ({percentage:.1f}%)\n")

        return self.earned_points, self.total_points, self.feedback


def main():
    """Main grading function."""
    if len(sys.argv) != 3:
        print("Usage: grade_assignment.py <assignment.md> <answer_key.json>")
        sys.exit(1)

    assignment_path = sys.argv[1]
    answer_key_path = sys.argv[2]

    grader = AssignmentGrader(assignment_path, answer_key_path)
    earned, total, feedback = grader.grade()

    # Print feedback
    for line in feedback:
        print(line)

    # Exit with success/failure code
    if earned >= total * 0.7:  # 70% passing
        sys.exit(0)
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
