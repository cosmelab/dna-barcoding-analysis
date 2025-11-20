# Lesson 1: Navigation - Finding Your Way Around

**Duration**: 30 minutes
**Topics**: `pwd`, `ls`, `cd`, `tree`

---

## Introduction

In this lesson, you'll learn to navigate the Linux filesystem. Think of it like learning to walk around a building - you need to know where you are, what's around you, and how to get to different rooms.

**Key Concept**: Linux filesystems are organized as a **tree structure**, starting from the root directory `/`.

---

## The Linux Filesystem

```
/                           # Root (top of everything)
├── home/                   # User home directories
│   ├── student1/
│   ├── student2/
│   └── luciano/
├── usr/                    # User programs
├── etc/                    # Configuration files
├── var/                    # Variable data (logs, etc.)
└── tmp/                    # Temporary files
```

**Important Directories**:
- `/` - Root directory (top of the tree)
- `/home/username/` - Your personal directory
- `~` - Shortcut for your home directory
- `.` - Current directory
- `..` - Parent directory (one level up)

---

## Command 1: `pwd` (Print Working Directory)

**Purpose**: Shows you where you are right now.

**Syntax**:
```bash
pwd
```

**Example**:
```bash
$ pwd
/home/luciano/dna-barcoding-analysis
```

**When to Use**:
- When you're lost and need to know your location
- Before running commands that affect your current directory
- In scripts to verify you're in the right place

**Try It**:
```bash
# Open terminal and type:
pwd

# You should see your current directory path
```

---

## Command 2: `ls` (List)

**Purpose**: Shows you what files and directories are in your current location.

**Basic Syntax**:
```bash
ls                  # List files in current directory
ls directory_name   # List files in specified directory
```

**Common Options**:
```bash
ls -l               # Long format (detailed info)
ls -a               # Show all files (including hidden)
ls -h               # Human-readable file sizes
ls -t               # Sort by modification time
ls -r               # Reverse order
ls -R               # Recursive (show subdirectories)

# Combine options
ls -lah             # Long format, all files, human-readable
ls -ltr             # Long format, sorted by time, reversed (oldest first)
```

**Examples**:

```bash
# Basic list
$ ls
00_introduction  01_linux_basics  02_python_basics  README.md

# Long format with details
$ ls -l
drwxr-xr-x  5 luciano  staff   160 Nov  7 14:30 00_introduction
drwxr-xr-x  8 luciano  staff   256 Nov  7 15:45 01_linux_basics
drwxr-xr-x  6 luciano  staff   192 Nov  7 14:30 02_python_basics
-rw-r--r--  1 luciano  staff  8432 Nov  7 16:20 README.md

# Show hidden files (start with .)
$ ls -a
.  ..  .git  .gitignore  .tracking  00_introduction  01_linux_basics

# Human-readable sizes
$ ls -lh
-rw-r--r--  1 luciano  staff   8.2K Nov  7 16:20 README.md

# Sort by modification time
$ ls -lt
# (shows most recently modified files first)
```

**Understanding `ls -l` Output**:
```
drwxr-xr-x  5  luciano  staff   160  Nov  7 14:30  00_introduction
│           │    │        │      │      │            │
│           │    │        │      │      │            └─ Filename
│           │    │        │      │      └─ Modification date/time
│           │    │        │      └─ File size (bytes)
│           │    │        └─ Group
│           │    └─ Owner
│           └─ Number of links
└─ Permissions (d=directory, r=read, w=write, x=execute)
```

**Try It**:
```bash
# List files in current directory
ls

# List with details
ls -l

# List all files including hidden
ls -a

# List with human-readable sizes
ls -lh

# List a specific directory
ls 01_linux_basics
```

---

## Command 3: `cd` (Change Directory)

**Purpose**: Move to a different directory.

**Syntax**:
```bash
cd directory_name       # Go to directory
cd ..                   # Go up one level
cd ../..                # Go up two levels
cd ~                    # Go to home directory
cd -                    # Go to previous directory
cd                      # Go to home (same as cd ~)
cd /                    # Go to root directory
```

**Absolute vs. Relative Paths**:

**Absolute Path** (starts with `/`):
```bash
cd /home/luciano/dna-barcoding-analysis/01_linux_basics
# Goes to exact location regardless of where you are
```

**Relative Path** (relative to current location):
```bash
cd 01_linux_basics      # Go to subdirectory
cd ../02_python_basics  # Go up one level, then into sibling directory
```

**Examples**:

```bash
# Start in home directory
$ pwd
/home/luciano

# Go to a subdirectory
$ cd dna-barcoding-analysis
$ pwd
/home/luciano/dna-barcoding-analysis

# Go to a subdirectory of current directory
$ cd 01_linux_basics
$ pwd
/home/luciano/dna-barcoding-analysis/01_linux_basics

# Go up one level
$ cd ..
$ pwd
/home/luciano/dna-barcoding-analysis

# Go up two levels
$ cd ../..
$ pwd
/home/luciano

# Go to home directory (two ways)
$ cd ~
# or just
$ cd

# Go back to previous directory
$ cd -
/home/luciano/dna-barcoding-analysis

# Go to root
$ cd /
$ pwd
/
```

**Try It**:
```bash
# Save your starting location
pwd

# Go to home directory
cd ~
pwd

# Go to root
cd /
ls

# Go back home
cd ~

# Navigate to the dna-barcoding-analysis directory
cd dna-barcoding-analysis
ls

# Go into 01_linux_basics
cd 01_linux_basics
pwd

# Go back up
cd ..
pwd

# Go back to where you started
cd -
```

---

## Command 4: `tree` (Display Directory Tree)

**Purpose**: Shows directory structure as a tree (if installed).

**Syntax**:
```bash
tree                    # Show tree of current directory
tree -L 2               # Limit depth to 2 levels
tree -d                 # Directories only
tree directory_name     # Show tree of specific directory
```

**Example**:
```bash
$ tree -L 2 .
.
├── 00_introduction
│   └── README.md
├── 01_linux_basics
│   ├── README.md
│   ├── 01_navigation.md
│   ├── data
│   ├── exercises
│   └── solutions
├── 02_python_basics
│   ├── README.md
│   └── notebooks
└── README.md
```

**Note**: If `tree` is not installed, use:
```bash
# Alternative: use ls with recursion
ls -R

# Or find (covered later)
find . -type d
```

---

## Practical Examples for Bioinformatics

### Example 1: Navigate to Your Sequences
```bash
# Where am I?
pwd

# Go to data directory
cd 04_data/student_sequences

# List my sequence files
ls -lh

# Go back to project root
cd ../../
```

### Example 2: Explore Module Structure
```bash
# See all modules
ls -d */

# Go into alignment module
cd 06_alignment

# What's in here?
ls -l

# Check examples subdirectory
ls examples/

# Go back
cd ..
```

### Example 3: Find Large Files
```bash
# List files sorted by size
ls -lhS

# Show sizes in human-readable format
ls -lh
```

---

## Tips and Tricks

### 1. Tab Completion is Your Friend
```bash
# Type first few letters, then press Tab
cd 01_li[TAB]
# Expands to: cd 01_linux_basics/
```

### 2. Use Up/Down Arrows
- Press Up Arrow to see previous commands
- Press Down Arrow to go forward in history
- Saves retyping commands

### 3. Use `cd -` to Toggle
```bash
cd /long/path/to/directory1
# Do some work
cd /different/long/path/to/directory2
# Do some work
cd -  # Back to directory1
cd -  # Back to directory2
```

### 4. Wildcards with `ls`
```bash
ls *.fasta              # List all FASTA files
ls sample_*             # List files starting with "sample_"
ls -d */                # List only directories
```

### 5. Check Before You Leap
```bash
# Before cd, use ls to see what's there
ls target_directory
cd target_directory
```

---

## Common Errors and Solutions

### Error 1: "No such file or directory"
```bash
$ cd nonexistent_folder
bash: cd: nonexistent_folder: No such file or directory

# Solution: Check spelling, use ls to verify
ls
```

### Error 2: "Permission denied"
```bash
$ cd /root
bash: cd: /root: Permission denied

# Solution: You don't have access to that directory
# Use a directory you have access to
```

### Error 3: Spaces in Directory Names
```bash
# WRONG
$ cd my folder
bash: cd: my: No such file or directory

# RIGHT (use quotes)
$ cd "my folder"

# BETTER (avoid spaces, use underscores)
$ cd my_folder
```

---

## Practice Exercises

### Exercise 1: Basic Navigation

```bash
# 1. Print your current directory
pwd

# 2. Go to your home directory
cd ~

# 3. List all files including hidden ones
ls -a

# 4. Go to the root directory
cd /

# 5. List only directories
ls -d */

# 6. Go back to home
cd ~

# 7. Navigate to dna-barcoding-analysis
cd dna-barcoding-analysis

# 8. List all modules
ls -d */

# 9. Go into 01_linux_basics
cd 01_linux_basics

# 10. Go back up one level
cd ..
```

### Exercise 2: Relative vs. Absolute Paths

```bash
# Starting from dna-barcoding-analysis directory

# 1. Use relative path to go to 02_python_basics
cd 02_python_basics

# 2. Use relative path to go to 03_r_basics
cd ../03_r_basics

# 3. Use absolute path to go to 01_linux_basics
cd /full/path/to/dna-barcoding-analysis/01_linux_basics

# 4. Use cd - to go back
cd -

# 5. Go to home using shortcut
cd ~
```

### Exercise 3: Exploring the Repository

```bash
# 1. Navigate to dna-barcoding-analysis
cd ~/dna-barcoding-analysis

# 2. How many top-level directories are there?
ls -d */ | wc -l

# 3. List all README files
ls */README.md

# 4. Go to the data directory
cd 04_data

# 5. What subdirectories exist?
ls -d */

# 6. Go back to project root
cd ..

# 7. Show the directory tree (2 levels deep)
tree -L 2
# or
ls -R | head -30
```

---

## Quick Reference

```bash
pwd                     # Where am I?
ls                      # What's here?
ls -lah                 # Show all files with details
cd directory            # Go somewhere
cd ..                   # Go up one level
cd ~                    # Go home
cd -                    # Go to previous directory
tree -L 2               # Show directory tree
```

---

## Check Your Understanding

Can you answer these without looking?

1. What command shows your current directory?
2. What does `ls -a` do?
3. How do you go up one directory level?
4. What's the shortcut for your home directory?
5. How do you go to the previous directory?
6. What's the difference between `cd /data` and `cd data`?

<details>
<summary>Click to see answers</summary>

1. `pwd`
2. Shows all files, including hidden files (starting with .)
3. `cd ..`
4. `~` (tilde)
5. `cd -`
6. `/data` is absolute path (starts from root), `data` is relative path (from current location)

</details>

---

## What's Next?

Now that you can navigate the filesystem, proceed to:

**Lesson 2: File Manipulation** - Learn to create, copy, move, and delete files.

**File**: `02_file_manipulation.md`

---

**Navigation mastered? Move on to Lesson 2!** →
