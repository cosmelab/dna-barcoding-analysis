#!/usr/bin/env python3
"""
Data Loader - Lab Data Analysis Module

Loads all class data CSVs and provides unified access.
Uses student codes only (HV, JR, MA, etc.) - no names or netids.
"""

import pandas as pd
from pathlib import Path


class LabData:
    """Loads and provides access to all lab data."""

    def __init__(self, data_dir=None):
        """Initialize with data directory path."""
        if data_dir is None:
            # Auto-detect relative to this file
            script_dir = Path(__file__).parent
            workspace_root = script_dir.parent.parent
            self.data_dir = workspace_root / "data" / "class_data"
        else:
            self.data_dir = Path(data_dir)

        # Source data directory (for gel images)
        self.source_dir = self.data_dir.parent.parent / "student_data_raw" / "source_data"

        # Data storage
        self.raw_measurements = None
        self.column_extraction = None
        self.pcr_summary = None
        self.nanodrop = None
        self.quality_assessment = None
        self.gel_results = None
        self.wolbachia = None
        self.sequencing = None
        self.teams = None

    def load_all(self):
        """Load all data files."""
        print("Loading lab data...")

        # HMW DNA extraction measurements
        f = self.data_dir / "raw_measurements.csv"
        if f.exists():
            self.raw_measurements = pd.read_csv(f)
            self._standardize_student_col(self.raw_measurements)
            print(f"  ✓ raw_measurements: {len(self.raw_measurements)} rows")

        # Column extraction (Zymo kit)
        f = self.data_dir / "column_extraction.csv"
        if f.exists():
            self.column_extraction = pd.read_csv(f)
            self._standardize_student_col(self.column_extraction)
            print(f"  ✓ column_extraction: {len(self.column_extraction)} rows")

        # PCR products summary
        f = self.data_dir / "pcr_products_summary.csv"
        if f.exists():
            self.pcr_summary = pd.read_csv(f)
            self._standardize_student_col(self.pcr_summary)
            print(f"  ✓ pcr_summary: {len(self.pcr_summary)} rows")

        # NanoDrop quality ratios
        f = self.data_dir / "nanodrop_quality.csv"
        if f.exists():
            self.nanodrop = pd.read_csv(f)
            self._standardize_student_col(self.nanodrop)
            print(f"  ✓ nanodrop: {len(self.nanodrop)} rows")

        # Quality assessment grades
        f = self.data_dir / "extraction_quality_assessment.csv"
        if f.exists():
            self.quality_assessment = pd.read_csv(f)
            self._standardize_student_col(self.quality_assessment)
            print(f"  ✓ quality_assessment: {len(self.quality_assessment)} rows")

        # COI PCR gel results
        f = self.data_dir / "coi_pcr_gel_results.csv"
        if f.exists():
            self.gel_results = pd.read_csv(f)
            self._standardize_student_col(self.gel_results)
            print(f"  ✓ gel_results: {len(self.gel_results)} rows")

        # Wolbachia PCR
        f = self.data_dir / "wolbachia_pcr.csv"
        if f.exists():
            self.wolbachia = pd.read_csv(f)
            self._standardize_student_col(self.wolbachia)
            print(f"  ✓ wolbachia: {len(self.wolbachia)} rows")

        # Sequencing results
        f = self.data_dir / "sequencing_success_summary.csv"
        if f.exists():
            self.sequencing = pd.read_csv(f)
            self._standardize_student_col(self.sequencing)
            print(f"  ✓ sequencing: {len(self.sequencing)} rows")

        # Team assignments
        f = self.source_dir / "teams.csv"
        if f.exists():
            self.teams = pd.read_csv(f, usecols=['Student', 'Team'])
            self._standardize_student_col(self.teams)
            print(f"  ✓ teams: {len(self.teams)} rows")

        print("Done.\n")
        return self

    def _standardize_student_col(self, df):
        """Rename 'Student' to 'Student_Code' if present."""
        if 'Student' in df.columns and 'Student_Code' not in df.columns:
            df.rename(columns={'Student': 'Student_Code'}, inplace=True)

    def get_gel_images(self):
        """Return paths to gel image files."""
        gels = {
            'coi_gel1': self.source_dir / "gel1_PCR_COI_WL-KG-MA-HV.png",
            'coi_gel2': self.source_dir / "gel2_PCR_COI_TW-BR-WA-JR.png",
            'coi_gel3': self.source_dir / "gel3_PCR_COI_JA-JM.png",
            'wolbachia': self.source_dir / "wolbachia_gel.png"
        }
        return {k: v for k, v in gels.items() if v.exists()}

    def get_team(self, student_code):
        """Get team assignment for a student."""
        if self.teams is None:
            return None
        row = self.teams[self.teams['Student_Code'] == student_code]
        if len(row) > 0:
            return row.iloc[0]['Team']
        return None

    def get_students_by_team(self, team):
        """Get list of students in a team."""
        if self.teams is None:
            return []
        return self.teams[self.teams['Team'] == team]['Student_Code'].tolist()


if __name__ == "__main__":
    # Test loading
    data = LabData().load_all()

    print("Gel images found:")
    for name, path in data.get_gel_images().items():
        print(f"  {name}: {path.name}")

    print("\nTeams:")
    print(f"  Spin: {data.get_students_by_team('Spin')}")
    print(f"  Magnet: {data.get_students_by_team('Magnet')}")
