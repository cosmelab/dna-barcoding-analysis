#!/usr/bin/env python3
"""
Trim reference sequences to COI barcode region (~712bp)
Removes sequences outside the AUCOS amplicon region
"""

from Bio import SeqIO
from Bio.Seq import Seq
import sys
from pathlib import Path

# AUCOS primers from Hoque et al 2022 (degenerate bases handled)
FORWARD_PRIMER = "TATTTTCWACAAATCATAARGATATTGGWAC"  # W=A/T, R=A/G
REVERSE_PRIMER = "TAWACTTCWGGRTGWCCRAARAATCA"  # Reverse complement needed

def degenerate_match(primer, seq, max_mismatches=3):
    """Find primer in sequence allowing for degenerate bases"""
    degenerate = {
        'W': 'AT', 'R': 'AG', 'Y': 'CT', 'S': 'CG',
        'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT',
        'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'
    }

    primer_len = len(primer)
    for i in range(len(seq) - primer_len + 1):
        mismatches = 0
        for j, p_base in enumerate(primer):
            s_base = seq[i + j]
            if p_base in degenerate:
                if s_base not in degenerate[p_base]:
                    mismatches += 1
            elif p_base != s_base:
                mismatches += 1

            if mismatches > max_mismatches:
                break

        if mismatches <= max_mismatches:
            return i

    return -1

def find_barcode_region(sequence):
    """Find and extract the COI barcode region"""
    seq_str = str(sequence.upper())

    # Try to find forward primer
    fwd_pos = degenerate_match(FORWARD_PRIMER, seq_str)

    # Try reverse complement for reverse primer
    rev_comp = str(Seq(REVERSE_PRIMER).reverse_complement())
    rev_pos = degenerate_match(rev_comp, seq_str)

    if fwd_pos >= 0 and rev_pos >= 0 and rev_pos > fwd_pos:
        # Extract region between primers (inclusive)
        barcode = seq_str[fwd_pos:rev_pos + len(rev_comp)]
        return barcode, f"trimmed_{fwd_pos}_{rev_pos}"

    # If primers not found but sequence is reasonable length, keep it
    if 640 <= len(seq_str) <= 800:
        return seq_str, "kept_original"

    # If too long and can't find primers, trim to ~700bp from start
    # (many sequences start at barcode region)
    if len(seq_str) > 800:
        return seq_str[:750], "trimmed_front_750bp"

    # Too short
    return None, "too_short"

def main():
    if len(sys.argv) < 3:
        print("Usage: python trim_references_to_barcode.py <input.fasta> <output.fasta>")
        sys.exit(1)

    input_file = Path(sys.argv[1])
    output_file = Path(sys.argv[2])

    print(f"Trimming COI sequences to barcode region (~712bp)")
    print(f"Input: {input_file}")
    print(f"Output: {output_file}")
    print()

    trimmed_seqs = []
    stats = {"trimmed": 0, "kept": 0, "removed": 0}

    for record in SeqIO.parse(input_file, "fasta"):
        barcode, status = find_barcode_region(record.seq)

        if barcode:
            new_record = record[:]
            new_record.seq = Seq(barcode)
            trimmed_seqs.append(new_record)

            if "trimmed" in status:
                stats["trimmed"] += 1
                print(f"✂️  {record.id[:50]:50s} {len(record.seq):5d}bp -> {len(barcode):4d}bp")
            else:
                stats["kept"] += 1
                print(f"✓  {record.id[:50]:50s} {len(barcode):4d}bp (kept)")
        else:
            stats["removed"] += 1
            print(f"✗  {record.id[:50]:50s} {len(record.seq):5d}bp REMOVED ({status})")

    # Write output
    SeqIO.write(trimmed_seqs, output_file, "fasta")

    print()
    print("=" * 70)
    print(f"Total sequences processed: {stats['trimmed'] + stats['kept'] + stats['removed']}")
    print(f"  Trimmed to barcode region: {stats['trimmed']}")
    print(f"  Kept as-is (good length): {stats['kept']}")
    print(f"  Removed (too short): {stats['removed']}")
    print(f"  Output sequences: {len(trimmed_seqs)}")

    # Show length distribution
    lengths = [len(rec.seq) for rec in trimmed_seqs]
    if lengths:
        print()
        print(f"Length range: {min(lengths)} - {max(lengths)} bp")
        print(f"Average length: {sum(lengths)/len(lengths):.0f} bp")

if __name__ == "__main__":
    main()
