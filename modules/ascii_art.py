#!/usr/bin/env python3
"""
ASCII Art for DNA Barcoding Pipeline

Beautiful terminal visualizations for the analysis workflow.
"""

# Dracula color codes for terminal
COLORS = {
    'purple': '\033[38;5;141m',
    'cyan': '\033[38;5;117m',
    'green': '\033[38;5;84m',
    'orange': '\033[38;5;215m',
    'pink': '\033[38;5;212m',
    'yellow': '\033[38;5;228m',
    'red': '\033[38;5;203m',
    'white': '\033[38;5;255m',
    'gray': '\033[38;5;245m',
    'bold': '\033[1m',
    'reset': '\033[0m'
}


def colorize(text: str, color: str) -> str:
    """Add color to text."""
    return f"{COLORS.get(color, '')}{text}{COLORS['reset']}"


PIPELINE_BANNER = """
╔══════════════════════════════════════════════════════════════════════════════╗
║                                                                              ║
║      ██████╗ ███╗   ██╗ █████╗     ██████╗  █████╗ ██████╗  ██████╗ ██████╗  ║
║      ██╔══██╗████╗  ██║██╔══██╗    ██╔══██╗██╔══██╗██╔══██╗██╔════╝██╔═══██╗ ║
║      ██║  ██║██╔██╗ ██║███████║    ██████╔╝███████║██████╔╝██║     ██║   ██║ ║
║      ██║  ██║██║╚██╗██║██╔══██║    ██╔══██╗██╔══██║██╔══██╗██║     ██║   ██║ ║
║      ██████╔╝██║ ╚████║██║  ██║    ██████╔╝██║  ██║██║  ██║╚██████╗╚██████╔╝ ║
║      ╚═════╝ ╚═╝  ╚═══╝╚═╝  ╚═╝    ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝ ╚═════╝  ║
║                                                                              ║
║                     DNA Barcoding Analysis Pipeline                          ║
║                          ENTM 201L - UC Riverside                            ║
║                                                                              ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""


PIPELINE_WORKFLOW = """
┌──────────────────────────────────────────────────────────────────────────────┐
│                         ANALYSIS WORKFLOW                                    │
├──────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│   .ab1 FILES        STEP 1         STEP 2          STEP 3         STEP 4    │
│  ┌──────────┐    ┌──────────┐   ┌──────────┐   ┌──────────┐   ┌──────────┐  │
│  │ ████████ │    │    QC    │   │ CONSENSUS│   │  ALIGN   │   │   TREE   │  │
│  │ ████████ │───▶│  FILTER  │──▶│   F + R  │──▶│  MAFFT   │──▶│ IQ-TREE  │  │
│  │ ████████ │    │  Phred   │   │  merge   │   │          │   │          │  │
│  └──────────┘    └──────────┘   └──────────┘   └──────────┘   └──────────┘  │
│  Chromatograms   Quality Ctrl   Combine reads  Multi-align   Phylogenetics  │
│                                                                    │         │
│              ┌─────────────────────────────────────────────────────┘         │
│              │                                                               │
│              ▼                                                               │
│        ┌──────────┐   ┌──────────┐                                          │
│        │  BLAST   │   │   LAB    │                                          │
│        │   NCBI   │──▶│  REPORT  │                                          │
│        │ database │   │          │                                          │
│        └──────────┘   └──────────┘                                          │
│          STEP 5         STEP 6                                              │
│        Species ID    Visualizations                                         │
│                                                                              │
└──────────────────────────────────────────────────────────────────────────────┘
"""


PIPELINE_WORKFLOW_COMPACT = """
╭─────────────────────────────────────────────────────────────────────╮
│                    DNA BARCODING PIPELINE                           │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│  [.ab1]──▶[QC]──▶[CONSENSUS]──▶[ALIGN]──▶[TREE]──▶[BLAST]──▶[REPORT]│
│    │       │         │           │         │        │         │     │
│  Input   Filter    F+R pair    MAFFT    IQ-TREE   NCBI    Results   │
│                                                                     │
╰─────────────────────────────────────────────────────────────────────╯
"""


DNA_HELIX = """
        ╭───╮
       ╱ A-T ╲
      ╱       ╲
     │  G≡≡≡C  │
      ╲       ╱
       ╲ T-A ╱
        ╰───╯
       ╱ C≡≡≡G ╲
      ╱         ╲
     │   A-T    │
      ╲         ╱
       ╲ G≡≡≡C ╱
        ╰───╯
"""


MOSQUITO_ASCII = """
                           __
                          /  \\
            .---.        |    |
           /     \\  __   |    |
          |  o o  |/  \\  |    |
          |   <   |    \\ |    |
           \\_____/      \\|____|
            |   |    ~~~~~~~~~~~
           /|   |\\
          / |   | \\
         /  |   |  \\
            |   |
           /|   |\\
          '-'   '-'
"""


STEP_BOXES = {
    1: """
┌─────────────────────────────────────┐
│  STEP 1: QUALITY CONTROL            │
├─────────────────────────────────────┤
│                                     │
│   Input:  .ab1 chromatogram files   │
│   Check:  Phred scores (Q30+)       │
│   Check:  Sequence length (>500bp)  │
│   Output: passed_sequences.fasta    │
│                                     │
│   ▓▓▓▓▓▓▓▓▓▓▓░░░░  Quality Score    │
│   ACGTACGTACGT...  Good Sequence    │
│                                     │
└─────────────────────────────────────┘
""",
    2: """
┌─────────────────────────────────────┐
│  STEP 2: CONSENSUS SEQUENCES        │
├─────────────────────────────────────┤
│                                     │
│   Forward:  5'──ACGT...──────▶ 3'   │
│                    ║║║║              │
│   Reverse:  3'◀──────...TGCA──5'    │
│                    ║║║║              │
│   Consensus: ════ACGT...════        │
│                                     │
│   Combines F+R for 2x accuracy!     │
│                                     │
└─────────────────────────────────────┘
""",
    3: """
┌─────────────────────────────────────┐
│  STEP 3: SEQUENCE ALIGNMENT         │
├─────────────────────────────────────┤
│                                     │
│   Sample1:  ACGT--ACGTACGT          │
│   Sample2:  ACGTACACGTACGT          │
│   Sample3:  ACGT--ACGT-CGT          │
│   Ref:      ACGTACACGTACGT          │
│             ****  **** ***          │
│                                     │
│   MAFFT aligns all sequences        │
│                                     │
└─────────────────────────────────────┘
""",
    4: """
┌─────────────────────────────────────┐
│  STEP 4: PHYLOGENETIC TREE          │
├─────────────────────────────────────┤
│                                     │
│            ┌── Aedes aegypti        │
│         ┌──┤                        │
│         │  └── Aedes albopictus     │
│      ───┤                           │
│         │  ┌── Culex pipiens        │
│         └──┤                        │
│            └── YOUR SAMPLE ◄──      │
│                                     │
│   IQ-TREE builds ML phylogeny       │
│                                     │
└─────────────────────────────────────┘
""",
    5: """
┌─────────────────────────────────────┐
│  STEP 5: SPECIES IDENTIFICATION     │
├─────────────────────────────────────┤
│                                     │
│   Your sequence ──▶ NCBI GenBank    │
│                         │           │
│                         ▼           │
│   ┌─────────────────────────────┐   │
│   │ Top Hit: Culex pipiens      │   │
│   │ Identity: 99.2%             │   │
│   │ Accession: KY123456         │   │
│   └─────────────────────────────┘   │
│                                     │
│   BLAST searches global database    │
│                                     │
└─────────────────────────────────────┘
""",
    6: """
┌─────────────────────────────────────┐
│  STEP 6: LAB DATA ANALYSIS          │
├─────────────────────────────────────┤
│                                     │
│   DNA Extraction ──▶ PCR ──▶ Seq    │
│        │              │        │    │
│        ▼              ▼        ▼    │
│   ┌─────────┐   ┌────────┐  ┌────┐  │
│   │ Yield   │   │Success │  │ QC │  │
│   │ ▓▓▓▓░░░ │   │  60%   │  │40% │  │
│   └─────────┘   └────────┘  └────┘  │
│                                     │
│   Interactive Plotly visualizations │
│                                     │
└─────────────────────────────────────┘
"""
}


COMPLETION_BANNER = """
╔══════════════════════════════════════════════════════════════════════════════╗
║                                                                              ║
║                        ✓ ANALYSIS COMPLETE! ✓                                ║
║                                                                              ║
╠══════════════════════════════════════════════════════════════════════════════╣
║                                                                              ║
║   Your results are ready in the results/ directory:                          ║
║                                                                              ║
║   📊 01_qc/          Quality control report                                  ║
║   🧬 02_consensus/   Consensus sequences                                     ║
║   📐 03_alignment/   Sequence alignment                                      ║
║   🌳 04_phylogeny/   Phylogenetic tree                                       ║
║   🔍 05_blast/       Species identification                                  ║
║   📈 06_lab/         Lab data visualizations                                 ║
║                                                                              ║
║   Open the HTML reports in your browser to explore!                          ║
║                                                                              ║
╚══════════════════════════════════════════════════════════════════════════════╝
"""


ENVIRONMENT_CHOICE = """
╭─────────────────────────────────────────────────────────────────────╮
│                    CHOOSE YOUR ENVIRONMENT                          │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│   OPTION A: GitHub Codespaces (No Installation!)                    │
│   ──────────────────────────────────────────────                    │
│   • Click "Code" → "Codespaces" → "Create codespace"                │
│   • Everything runs in the cloud                                    │
│   • Use: ./tutorial-cs.sh and ./run-analysis-cs.sh                  │
│                                                                     │
│   OPTION B: Local Docker                                            │
│   ──────────────────────────────────────────────                    │
│   • Install Docker Desktop on your computer                         │
│   • Works on Mac (Intel/M1), Windows, Linux                         │
│   • Use: ./tutorial.sh and ./run-analysis.sh                        │
│                                                                     │
╰─────────────────────────────────────────────────────────────────────╯
"""


TEAM_BATTLE = """
╔══════════════════════════════════════════════════════════════════════════════╗
║                                                                              ║
║                    ⚔️  TEAM CHALLENGE: SPIN vs MAGNET  ⚔️                     ║
║                                                                              ║
╠═══════════════════════════════════╦══════════════════════════════════════════╣
║                                   ║                                          ║
║         TEAM SPIN                 ║           TEAM MAGNET                    ║
║       (Column-based)              ║         (Bead-based)                     ║
║                                   ║                                          ║
║    ┌───────────────────┐          ║      ┌───────────────────┐               ║
║    │   ◯◯◯◯◯◯◯◯◯◯     │          ║      │   ●●●●●●●●●●     │               ║
║    │   Spin columns    │          ║      │   Magnetic beads  │               ║
║    └───────────────────┘          ║      └───────────────────┘               ║
║                                   ║                                          ║
║    JR, HV, TW, JM, WL             ║      MA, BR, WA, KG, JA                  ║
║                                   ║                                          ║
╚═══════════════════════════════════╩══════════════════════════════════════════╝
"""


def print_banner():
    """Print the main pipeline banner."""
    print(colorize(PIPELINE_BANNER, 'purple'))


def print_workflow():
    """Print the workflow diagram."""
    print(colorize(PIPELINE_WORKFLOW, 'cyan'))


def print_workflow_compact():
    """Print compact workflow."""
    print(colorize(PIPELINE_WORKFLOW_COMPACT, 'cyan'))


def print_step(step_num: int):
    """Print a specific step box."""
    if step_num in STEP_BOXES:
        print(colorize(STEP_BOXES[step_num], 'green'))


def print_completion():
    """Print completion banner."""
    print(colorize(COMPLETION_BANNER, 'green'))


def print_environment_choice():
    """Print environment options."""
    print(colorize(ENVIRONMENT_CHOICE, 'yellow'))


def print_team_battle():
    """Print team challenge banner."""
    print(colorize(TEAM_BATTLE, 'orange'))


def print_dna_helix():
    """Print DNA helix."""
    print(colorize(DNA_HELIX, 'cyan'))


def print_mosquito():
    """Print mosquito ASCII art."""
    print(colorize(MOSQUITO_ASCII, 'gray'))


# Progress bar helper
def progress_bar(current: int, total: int, width: int = 40) -> str:
    """Create a progress bar string."""
    filled = int(width * current / total)
    bar = '▓' * filled + '░' * (width - filled)
    percent = current / total * 100
    return f"[{bar}] {percent:.0f}%"


if __name__ == "__main__":
    # Demo all ASCII art
    print_banner()
    print_workflow()
    print_environment_choice()

    for i in range(1, 7):
        print_step(i)

    print_team_battle()
    print_completion()

    print("\nProgress bar demo:")
    for i in range(0, 101, 20):
        print(f"  {progress_bar(i, 100)}")
