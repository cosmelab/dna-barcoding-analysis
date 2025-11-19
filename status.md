# Repository Status

**Last Updated**: 2025-11-18

---

## ✅ COMPLETED

### Working Analysis Modules (`modules/`)
- ✓ Module 00: F/R Assembly (script + README)
- ✓ Module 01: Quality Control (script + README)
- ✓ Module 02: Alignment (script + README)
- ✓ Module 03: Phylogeny (script + README)
- ✓ Module 04: Species ID (script + README)
- ✓ Master Pipeline (integrated all modules)

### Reference Data
- ✓ Hoque et al 2022 sequences fetched (52/53 mosquito COI)
- ✓ `fetch_reference_sequences.py` script

### Testing
- ✓ Tested with real UC genomics data (8 samples, 4 passed QC)
- ✓ F/R assembly working (AT_ROCK merged successfully)
- ✓ Chromatogram visualization working

### Documentation
- ✓ Module READMEs (00-04) - Student-friendly
- ✓ Tracking system (state.yaml, decisions.yaml, tasks.yaml, log.md)

---

## ⚠️ IN PROGRESS

### Docker Container
- ⚠️ Dockerfile created but **BUILD FAILING**
- Issue: colorls gem compilation (needs C compilers)
- Attempts: 7 builds, all failed
- Blockers: Missing dependencies for native gem extensions

### Tutorial Content (00-08 directories)
- ✓ READMEs exist (partially complete)
- ❌ Scripts subdirectories mostly empty
- ❌ Exercise solutions missing
- ❌ Example data incomplete

---

## ❌ NOT STARTED

### Tutorial Scripts (`scripts/`)
Placeholder directories exist but empty:
- `scripts/alignment/` - empty
- `scripts/phylogeny/` - empty
- `scripts/quality_control/` - empty
- `scripts/utilities/` - empty

### GitHub Classroom
- ❌ Template repository structure
- ❌ Autograding workflows
- ❌ Student starter files

### Additional Modules Needed
- ❌ Trimming module (user mentioned this)
- ❌ Wrapper scripts in /usr/local/bin for container

### Main Documentation
- ❌ Main README update
- ❌ Quick start guide
- ❌ Installation instructions
- ❌ Troubleshooting guide

---

## 🚨 CRITICAL BLOCKERS

### 1. Docker Build Failure
**Impact**: Students can't use the container
**Status**: 7 failed builds, last error: colorls gem native extension compilation
**Next Steps**:
- Option A: Keep debugging (add more dependencies)
- Option B: Remove colorls entirely (simpler terminal)
- Option C: Use pre-built colorls container layer

### 2. Scope Clarity
**Question**: Is this repo for:
- **Option A**: Week 8 lab only (just analysis modules)?
- **Option B**: Full quarter course (all 00-08 tutorials)?
- **Option C**: Both (comprehensive learning + analysis)?

---

## 📋 RECOMMENDED PRIORITY

### High Priority (Week 8 is soon)
1. **Fix Docker build** (critical blocker)
2. **Test full pipeline** with Hoque sequences
3. **Create wrapper scripts** for container commands
4. **Update main README** with quick start

### Medium Priority
5. **GitHub Classroom template**
6. **Populate scripts/** with examples
7. **Add trimming module**

### Low Priority (if time allows)
8. **Complete tutorial content** (00-08)
9. **Add exercise solutions**
10. **Create autograding workflows**

---

## 🎯 DECISION NEEDED

**What's the deadline and minimum viable product?**
- Is the container + modules/ enough for Week 8?
- Or do we need all tutorial content (00-08)?
- Can students work without container (local installation)?
