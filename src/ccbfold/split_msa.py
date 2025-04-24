from pathlib import Path
import click

@click.command()
@click.argument("filepath", type=Path)
def split_msa(filepath: Path) -> None:
    """Split a multiple sequence alignment file."""
    print(f"Splitting MSA hello at {filepath}")
