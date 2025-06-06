import click
from .wc import wc
from .split_msa import split_msa
from .eval_plots import plot_plddt, plot_pae

@click.group()
def main():
    """ccbfold CLI"""
    pass

main.add_command(wc)
main.add_command(split_msa)
