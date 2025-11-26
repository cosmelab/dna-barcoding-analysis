#!/usr/bin/env python3
"""
Automated Grading Script for DNA Barcoding Assignment

Reads answers.json and grades against answer_key.json
Since all students analyze the same class dataset, answers should be identical.
"""

import json
import sys
from pathlib import Path
from typing import Dict, List


class AnswerGrader:
    """Grade student answers from answers.json."""

    def __init__(self, answers_path: str, answer_key_path: str):
        self.answers_path = Path(answers_path)
        self.answer_key_path = Path(answer_key_path)
        self.earned_points = 0
        self.total_points = 40  # Part 2 (20) + Part 3 (20)
        self.feedback = []

    def load_json(self, path: Path) -> Dict:
        """Load JSON file."""
        with open(path) as f:
            return json.load(f)

    def grade_species_table(self, student: Dict, key: Dict) -> int:
        """Grade species identification table (20 points)."""
        points = 0
        max_points = 20

        student_table = student.get('species_table', [])
        key_table = key['species_table']

        if not student_table:
            self.feedback.append("❌ Results table is empty (0/20 points)")
            return 0

        if len(student_table) != len(key_table):
            self.feedback.append(
                f"⚠️  Expected {len(key_table)} samples, got {len(student_table)}"
            )

        # Grade each species entry
        points_per_sample = max_points / len(key_table)

        for i, expected in enumerate(key_table):
            if i >= len(student_table):
                self.feedback.append(f"❌ Missing entry for {expected['sample']}")
                continue

            student_entry = student_table[i]
            sample_points = 0

            # Check species name (case-insensitive)
            student_species = student_entry['species'].lower().strip()
            expected_species = expected['species'].lower().strip()

            if student_species == expected_species:
                sample_points += points_per_sample * 0.6  # 60% for species
                self.feedback.append(f"✓ {expected['sample']}: Correct species")
            else:
                self.feedback.append(
                    f"❌ {expected['sample']}: Wrong species. "
                    f"Expected '{expected['species']}', got '{student_entry['species']}'"
                )

            # Check % identity (allow ±1%)
            try:
                student_pct = float(str(student_entry['percent_identity']).replace('%', ''))
                expected_pct = float(str(expected['percent_identity']).replace('%', ''))

                if abs(student_pct - expected_pct) <= 1.0:
                    sample_points += points_per_sample * 0.4  # 40% for % identity
                    self.feedback.append(f"✓ {expected['sample']}: Correct % identity (~{student_pct}%)")
                else:
                    self.feedback.append(
                        f"❌ {expected['sample']}: Wrong % identity. "
                        f"Expected ~{expected_pct}%, got {student_pct}%"
                    )
            except (ValueError, KeyError):
                self.feedback.append(f"❌ {expected['sample']}: Invalid % identity format")

            points += sample_points

        self.feedback.append(f"\n✓ Results table: {points:.1f}/{max_points} points")
        return int(points)

    def grade_question1(self, student: Dict, key: Dict) -> int:
        """Grade Question 1 - Quality Control (5 points)."""
        points = 0
        max_points = 5

        # Q1a: Sequences passed (2 points)
        student_ans = str(student.get('q1a_sequences_passed', '')).strip()
        expected_ans = str(key['q1a_sequences_passed']).strip()

        if student_ans == expected_ans:
            points += 2
            self.feedback.append(f"✓ Q1(a): Correct ({expected_ans} sequences passed)")
        else:
            self.feedback.append(
                f"❌ Q1(a): Wrong answer. Expected {expected_ans}, got {student_ans}"
            )

        # Q1b: Pairs passed (2 points)
        student_ans = str(student.get('q1b_pairs_passed', '')).strip()
        expected_ans = str(key['q1b_pairs_passed']).strip()

        if student_ans == expected_ans:
            points += 2
            self.feedback.append(f"✓ Q1(b): Correct ({expected_ans} pairs)")
        else:
            self.feedback.append(
                f"❌ Q1(b): Wrong answer. Expected {expected_ans}, got {student_ans}"
            )

        # Q1c: Conceptual understanding (1 point)
        student_ans = student.get('q1c_explanation', '').lower()
        keywords = key['q1c_keywords']

        if any(kw in student_ans for kw in keywords):
            points += 1
            self.feedback.append("✓ Q1(c): Good explanation of F+R importance")
        else:
            self.feedback.append(
                f"❌ Q1(c): Missing key concepts. Should mention: {', '.join(keywords)}"
            )

        self.feedback.append(f"✓ Question 1: {points}/{max_points} points\n")
        return points

    def grade_question2(self, student: Dict, key: Dict) -> int:
        """Grade Question 2 - Phylogenetic Tree (7 points)."""
        points = 0
        max_points = 7

        # Q2a: Clustering (2 points)
        student_ans = student.get('q2a_clustering', '')
        expected_ans = key['q2a_clustering']

        if student_ans == expected_ans or any(kw in student_ans.lower() for kw in ['spread', 'different']):
            points += 2
            self.feedback.append("✓ Q2(a): Correct clustering description")
        else:
            self.feedback.append("❌ Q2(a): Incorrect clustering description")

        # Q2b: Related species (3 points)
        student_ans = student.get('q2b_related_species', '').lower()
        keywords = key['q2b_related_species_keywords']

        matches = sum(1 for kw in keywords if kw in student_ans)
        if matches >= len(keywords) * 0.5:  # At least 50% of keywords
            points += 3
            self.feedback.append("✓ Q2(b): Correctly identified related species")
        else:
            self.feedback.append(
                f"❌ Q2(b): Missing species. Should mention: {', '.join(keywords)}"
            )

        # Q2c: Diversity interpretation (2 points)
        student_ans = student.get('q2c_diversity', '').lower()
        keywords = key['q2c_diversity_keywords']

        if any(kw in student_ans for kw in keywords):
            points += 2
            self.feedback.append("✓ Q2(c): Good diversity interpretation")
        else:
            self.feedback.append("❌ Q2(c): Weak diversity interpretation")

        self.feedback.append(f"✓ Question 2: {points}/{max_points} points\n")
        return points

    def grade_question3(self, student: Dict, key: Dict) -> int:
        """Grade Question 3 - Species Identification (8 points)."""
        points = 0
        max_points = 8

        # Q3a: Species list (2 points)
        student_ans = student.get('q3a_species_list', '').lower()
        expected_species = [sp.lower() for sp in key['q3a_species_list']]

        matches = sum(1 for sp in expected_species if sp in student_ans)
        if matches >= len(expected_species) * 0.8:  # 80% of species
            points += 2
            self.feedback.append("✓ Q3(a): Correctly listed species")
        else:
            self.feedback.append(
                f"❌ Q3(a): Missing species. Expected: {', '.join(key['q3a_species_list'])}"
            )

        # Q3b: BLAST-tree agreement (2 points)
        student_ans = student.get('q3b_blast_tree_agree', '').lower()
        expected_ans = key['q3b_blast_tree_agree']

        if student_ans == expected_ans:
            points += 2
            self.feedback.append("✓ Q3(b): Correct BLAST-tree comparison")
        else:
            self.feedback.append("❌ Q3(b): Incorrect BLAST-tree comparison")

        # Q3c: SoCal occurrence (2 points)
        student_ans = student.get('q3c_socal_species', '').lower()
        expected_ans = key['q3c_socal_species']

        if student_ans == expected_ans:
            points += 2
            self.feedback.append("✓ Q3(c): Correct geographic assessment")
        else:
            self.feedback.append("❌ Q3(c): Incorrect geographic assessment")

        # Q3d: Confidence (2 points)
        student_ans = student.get('q3d_confidence', '')
        acceptable = key['q3d_confidence_acceptable']

        if any(acc in student_ans for acc in acceptable):
            points += 2
            self.feedback.append("✓ Q3(d): Appropriate confidence level")
        else:
            self.feedback.append("⚠️  Q3(d): Confidence level should be 'Very confident' or 'Confident'")

        self.feedback.append(f"✓ Question 3: {points}/{max_points} points\n")
        return points

    def grade(self) -> int:
        """Grade all answers and return score."""
        try:
            student_answers = self.load_json(self.answers_path)
            answer_key = self.load_json(self.answer_key_path)
        except FileNotFoundError as e:
            print(f"❌ Error: {e}")
            return 1

        self.feedback.append("\n" + "=" * 70)
        self.feedback.append("GRADING RESULTS")
        self.feedback.append("=" * 70 + "\n")

        # Grade each part
        self.earned_points += self.grade_species_table(student_answers, answer_key)
        self.earned_points += self.grade_question1(student_answers, answer_key)
        self.earned_points += self.grade_question2(student_answers, answer_key)
        self.earned_points += self.grade_question3(student_answers, answer_key)

        # Final score
        percentage = (self.earned_points / self.total_points) * 100

        self.feedback.append("=" * 70)
        self.feedback.append(f"FINAL SCORE: {self.earned_points}/{self.total_points} ({percentage:.1f}%)")
        self.feedback.append("=" * 70 + "\n")

        # Print all feedback
        for line in self.feedback:
            print(line)

        # Return exit code (0 = pass, 1 = fail)
        return 0 if percentage >= 70 else 1


def main():
    """Main entry point."""
    if len(sys.argv) != 3:
        print("Usage: grade_answers.py <answers.json> <answer_key.json>")
        sys.exit(1)

    grader = AnswerGrader(sys.argv[1], sys.argv[2])
    exit_code = grader.grade()
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
