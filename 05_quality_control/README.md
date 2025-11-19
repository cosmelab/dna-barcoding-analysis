# Module 05: Quality Control of Sanger Sequences

**Duration**: 2 hours
**Prerequisites**: Modules 01-02

---

## Overview

Learn to assess and clean DNA sequencing data from Sanger sequencing. Quality control is critical for reliable downstream analysis.

---

## Contents

```
05_quality_control/
├── README.md
├── scripts/
│   ├── assess_chromatograms.py      # Parse .ab1 files
│   ├── trim_sequences.sh            # Trim low quality ends
│   ├── quality_report.py            # Generate QC report
│   └── filter_sequences.py          # Filter by quality
└── examples/
    ├── good_quality.ab1
    ├── poor_quality.ab1
    └── contamination_example.ab1
```

---

## Quality Metrics

### Phred Quality Scores

```
Q = -10 × log₁₀(P)
```

- **Q30** = 99.9% accuracy (1 error per 1000 bases)
- **Q20** = 99% accuracy (1 error per 100 bases)
- **Q10** = 90% accuracy (1 error per 10 bases)

**Rule of Thumb**: Keep bases with Q > 20

---

## Workflow

### 1. Inspect Chromatograms

```python
from Bio import SeqIO

# Read .ab1 file
record = SeqIO.read("sample.ab1", "abi")

# Get quality scores
qualities = record.letter_annotations["phred_quality"]

# Check sequence
print(f"Sequence length: {len(record.seq)}")
print(f"Average quality: {sum(qualities)/len(qualities):.1f}")
```

### 2. Trim Low-Quality Ends

```bash
# Using seqtk
seqtk trimfq -q 0.01 input.fastq > trimmed.fastq

# Using BioPython
python scripts/trim_sequences.sh input_dir/ output_dir/
```

### 3. Check for Issues

**Common Problems**:
- Low quality at ends (normal, trim it)
- Double peaks (contamination or heterozygosity)
- Short reads (<500 bp after trimming)
- High N content (ambiguous bases)

---

## Example: Complete QC Pipeline

```python
#!/usr/bin/env python3
"""
Quality control pipeline for Sanger sequences
"""

from Bio import SeqIO
import sys
import os

def assess_quality(ab1_file, min_quality=20, min_length=500):
    """Assess sequence quality."""
    record = SeqIO.read(ab1_file, "abi")
    qualities = record.letter_annotations["phred_quality"]

    # Calculate metrics
    avg_quality = sum(qualities) / len(qualities)
    bases_above_q20 = sum(1 for q in qualities if q >= 20)
    percent_q20 = (bases_above_q20 / len(qualities)) * 100

    # Pass/Fail
    passed = (avg_quality >= min_quality and
              len(record.seq) >= min_length and
              percent_q20 >= 90)

    return {
        "file": os.path.basename(ab1_file),
        "length": len(record.seq),
        "avg_quality": avg_quality,
        "percent_q20": percent_q20,
        "passed": passed
    }

def main():
    import glob

    # Process all .ab1 files
    results = []
    for ab1_file in glob.glob("*.ab1"):
        result = assess_quality(ab1_file)
        results.append(result)
        print(f"{result['file']}: {'PASS' if result['passed'] else 'FAIL'}")

    # Summary
    passed = sum(1 for r in results if r["passed"])
    print(f"\nSummary: {passed}/{len(results)} sequences passed QC")

if __name__ == "__main__":
    main()
```

---

## Quality Flags

| Flag | Meaning | Action |
|------|---------|--------|
| ✓ | Good quality (Q>30) | Proceed to alignment |
| ⚠ | Moderate quality (Q20-30) | Use with caution |
| ✗ | Poor quality (Q<20) | Re-sequence or discard |
| N | Ambiguous base | Remove or trim |
| ? | Double peak | Check for contamination |

---

## Tools

- **BioPython**: Parse .ab1 files
- **seqtk**: Fast trimming
- **Phred/Phrap**: Classic QC tools
- **FastQC**: Quality visualization (for NGS, but works for Sanger too)

---

## Next Steps

Clean sequences proceed to **Module 06: Alignment** →
