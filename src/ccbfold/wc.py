from pathlib import Path
import click

def wc_impl(filepath: Path, show_bytes: bool, show_words: bool, show_lines: bool) -> None:
    content = filepath.read_bytes()
    byte_count = len(content)
    word_count = len(content.split())
    line_count = len(content.splitlines())
    output = (
        (f"{line_count:8}" if show_lines else "")
        + (f"{word_count:8}" if show_words else "")
        + (f"{byte_count:8}" if show_bytes else "")
        + f" {filepath}"
    )
    print(output)

@click.command()
@click.argument("filepath", type=Path)
@click.option("-c", is_flag=True, help="Show byte count.")
@click.option("-w", is_flag=True, help="Show word count.")
@click.option("-l", is_flag=True, help="Show line count.")
def wc(filepath: Path, c: bool, w: bool, l: bool) -> None:
    """Count lines, words, and bytes in a file."""
    if {c, w, l} == {False}:
        c, w, l = True, True, True
    wc_impl(filepath, c, w, l)
