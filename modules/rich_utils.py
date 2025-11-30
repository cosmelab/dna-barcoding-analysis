#!/usr/bin/env python3
"""
Rich Terminal Utilities for DNA Barcoding Pipeline

Provides beautiful terminal output with Rich library.
Falls back to plain text if Rich is not available.
"""

try:
    from rich.console import Console
    from rich.panel import Panel
    from rich.table import Table
    from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn
    from rich.status import Status
    from rich.text import Text
    from rich import box
    RICH_AVAILABLE = True
except ImportError:
    RICH_AVAILABLE = False
    Status = None

# Global console instance with force_terminal=True
# This ensures colors work even when output is piped through tee
console = Console(force_terminal=True) if RICH_AVAILABLE else None


def print_header(title: str, subtitle: str = ""):
    """Print a beautiful header panel."""
    if RICH_AVAILABLE:
        content = f"[bold cyan]{title}[/]"
        if subtitle:
            content += f"\n[dim]{subtitle}[/]"
        console.print()
        console.print(Panel(content, box=box.DOUBLE, border_style="cyan"))
        console.print()
    else:
        print("\n" + "=" * 70)
        print(f"  {title}")
        if subtitle:
            print(f"  {subtitle}")
        print("=" * 70 + "\n")


def print_step(step_num: int, total_steps: int, description: str):
    """Print a step indicator."""
    if RICH_AVAILABLE:
        console.print()
        console.rule(f"[bold yellow]Step {step_num}/{total_steps}: {description}[/]")
        console.print()
    else:
        print(f"\nStep {step_num}/{total_steps}: {description}")
        print("-" * 70)


def print_success(message: str):
    """Print a success message."""
    if RICH_AVAILABLE:
        console.print(f"[bold green]✓[/] {message}")
    else:
        print(f"✓ {message}")


def print_error(message: str):
    """Print an error message."""
    if RICH_AVAILABLE:
        console.print(f"[bold red]✗[/] {message}")
    else:
        print(f"✗ {message}")


def print_warning(message: str):
    """Print a warning message."""
    if RICH_AVAILABLE:
        console.print(f"[yellow]⚠[/] {message}")
    else:
        print(f"⚠ {message}")


def print_info(message: str, indent: bool = True):
    """Print an info message."""
    prefix = "  " if indent else ""
    if RICH_AVAILABLE:
        console.print(f"{prefix}[dim]{message}[/]")
    else:
        print(f"{prefix}{message}")


def print_file(label: str, path: str):
    """Print a file path nicely."""
    if RICH_AVAILABLE:
        console.print(f"[green]✓[/] {label}: [cyan]{path}[/]")
    else:
        print(f"✓ {label}: {path}")


def print_summary(title: str, items: dict):
    """Print a summary table."""
    if RICH_AVAILABLE:
        table = Table(title=title, box=box.ROUNDED, show_header=False)
        table.add_column("Metric", style="cyan")
        table.add_column("Value", style="bold")

        for key, value in items.items():
            table.add_row(key, str(value))

        console.print()
        console.print(table)
        console.print()
    else:
        print(f"\n{title}")
        print("-" * 40)
        for key, value in items.items():
            print(f"  {key}: {value}")
        print()


def print_complete(module_name: str, output_path: str = ""):
    """Print completion message."""
    if RICH_AVAILABLE:
        content = f"[bold green]{module_name} Complete![/]"
        if output_path:
            content += f"\n\n[dim]Output:[/] [cyan]{output_path}[/]"
        console.print()
        console.print(Panel(content, border_style="green"))
        console.print()
    else:
        print(f"\n{'=' * 70}")
        print(f"  {module_name} COMPLETE")
        print("=" * 70)
        if output_path:
            print(f"Output: {output_path}")
        print()


def create_progress():
    """Create a progress bar context manager."""
    if RICH_AVAILABLE:
        return Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            console=console
        )
    else:
        return None


class SimpleProgress:
    """Simple progress indicator for non-Rich environments."""
    def __init__(self, total: int, description: str = "Processing"):
        self.total = total
        self.current = 0
        self.description = description

    def update(self, n: int = 1):
        self.current += n
        pct = (self.current / self.total) * 100
        print(f"\r  {self.description}: {self.current}/{self.total} ({pct:.0f}%)", end="", flush=True)
        if self.current >= self.total:
            print()  # Newline at end


def get_spinner(message: str):
    """Get a spinner context manager for long-running operations."""
    if RICH_AVAILABLE and console:
        return console.status(f"[bold cyan]{message}[/]", spinner="dots")
    else:
        # Return a dummy context manager for non-Rich
        class DummySpinner:
            def __enter__(self):
                print(f"  {message}...")
                return self
            def __exit__(self, *args):
                pass
            def update(self, msg):
                pass
        return DummySpinner()
