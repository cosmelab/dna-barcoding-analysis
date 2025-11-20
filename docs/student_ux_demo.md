# Student User Experience - Terminal Output Demo

This document shows what students will see when running the analysis pipeline.

## Improved Features

1. **Clear step-by-step progress** - Students know exactly what's happening
2. **Visual feedback** - Checkmarks, progress indicators, section headers
3. **Helpful error messages** - Clear guidance when something goes wrong
4. **File locations displayed** - Easy to find results
5. **Next steps guidance** - Students know what to do after each step
6. **Auto-open HTML** - Optional `--open` flag to view results immediately

## Example: Quality Control Module

### Command
```bash
python qc_chromatograms.py data/my_sequences/ results/ --open
```

### Terminal Output
```
======================================================================
  DNA BARCODING QUALITY CONTROL
======================================================================
Analyzing Sanger chromatograms (.ab1 files)
Input: data/my_sequences
Output: results

[Step 1/3] Finding chromatogram files
----------------------------------------------------------------------
✓ Found 8 chromatogram files
  • AT83F_A01_015.ab1
  • AT83R_E01_007.ab1
  • AT94F_B01_013.ab1
  • AT94R_F01_005.ab1
  • AT99F_C01_011.ab1
  • AT99R_G01_003.ab1
  • AT_ROCK_F_D01_009.ab1
  • AT_ROCK_R_H01_001.ab1

[Step 2/3] Analyzing chromatograms
----------------------------------------------------------------------
  Checking sequence quality, length, and reading frames...
  [1/8] AT83F_A01_015.ab1...
      ✗ FAILED (Low quality - avg Q score: 18.5)
  [2/8] AT83R_E01_007.ab1...
      ✗ FAILED (Too short - 412 bp)
  [3/8] AT94F_B01_013.ab1...
      ✗ FAILED (Low quality - avg Q score: 19.2)
  [4/8] AT94R_F01_005.ab1...
      ✗ FAILED (Low quality - avg Q score: 17.8)
  [5/8] AT99F_C01_011.ab1...
      ✓ PASSED (length: 687 bp, avg quality: 47.3)
  [6/8] AT99R_G01_003.ab1...
      ✓ PASSED (length: 672 bp, avg quality: 52.5)
  [7/8] AT_ROCK_F_D01_009.ab1...
      ✓ PASSED (length: 636 bp, avg quality: 44.0)
  [8/8] AT_ROCK_R_H01_001.ab1...
      ✓ PASSED (length: 688 bp, avg quality: 42.8)

[Step 3/3] Generating reports
----------------------------------------------------------------------
  Creating HTML report with chromatogram visualizations...
✓ HTML report: results/qc_report.html
✓ JSON data: results/qc_results.json
✓ Passed sequences: results/passed_sequences.fasta

======================================================================
  QUALITY CONTROL COMPLETE
======================================================================
Total sequences analyzed: 8
✓ Passed QC: 4 sequences
✗ Failed QC: 4 sequences

======================================================================
  NEXT STEPS:
======================================================================
1. Open the HTML report to view detailed results:
   /Users/student/dna-barcoding/results/qc_report.html

2. Passed sequences are ready for alignment:
   /Users/student/dna-barcoding/results/passed_sequences.fasta
======================================================================

Opening HTML report in your web browser...
✓ Report opened successfully
```

## Cross-Platform Browser Opening

The `--open` flag uses Python's `webbrowser` module which works on:
- **Mac** - Opens in default browser (Safari, Chrome, Firefox, etc.)
- **Windows** - Opens in default browser (Edge, Chrome, Firefox, etc.)
- **Linux** - Opens in default browser (Firefox, Chrome, etc.)

### How it works:
1. Converts relative paths to absolute paths
2. Creates proper `file://` URLs
3. Calls system's default browser
4. Falls back gracefully if browser can't be opened

## Benefits for Students

### Before (old output):
```
Found 8 chromatogram files
Analyzing AT83F_A01_015.ab1...
Analyzing AT83R_E01_007.ab1...
...
HTML report: /workspace/results/qc_report.html/qc_report.html
Summary:
  Total: 8
  Passed: 4
```

Students didn't know:
- What step they were on
- Why sequences failed
- Where to find results
- What to do next

### After (new output):
```
[Step 2/3] Analyzing chromatograms
  [5/8] AT99F_C01_011.ab1...
      ✓ PASSED (length: 687 bp, avg quality: 47.3)
```

Students can see:
- ✓ Current step (2 of 3)
- ✓ Progress (5 of 8 files)
- ✓ What's happening (analyzing chromatograms)
- ✓ Results (passed/failed with reasons)
- ✓ Quality metrics (length, quality score)
- ✓ Next steps (clear instructions)
- ✓ Can auto-open results with `--open`

## Planned Improvements for All Modules

All modules will be updated with:
1. Consistent terminal output formatting
2. Step-by-step progress indicators
3. `--open` flag for viewing HTML reports
4. Clear "NEXT STEPS" guidance
5. Helpful error messages
6. Cross-platform compatibility

Modules to update:
- [x] `01_quality_control/qc_chromatograms.py` - DONE
- [ ] `02_alignment/align_sequences.py`
- [ ] `03_phylogeny/build_tree.py`
- [ ] `04_identification/identify_species.py`
- [ ] `master_pipeline.py`
