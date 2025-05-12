import json
from collections import Counter
from pathlib import Path

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np


def plot_plddt(atom_plddts: list[int], atom_chain_ids: list[str]) -> None:
    ''' Plot pLDDT data '''

    xlim = len(atom_plddts)+1

    fig = plt.figure(figsize=(15, 5))
    ax = fig.gca()

    ax.add_patch(
        patches.Rectangle(
            (1, 90), xlim, 10,
            color='#0751d9', alpha=1, label='Very High (90-100)'
        )
    )
    ax.add_patch(
        patches.Rectangle(
            (1, 70), xlim, 20,
            color='#69c9f4', alpha=1, label='High (70-90)'
        )
    )
    ax.add_patch(
        patches.Rectangle(
            (1, 50), xlim, 20,
            color='#fedb0f', alpha=1, label='Low (50-70)'
        )
    )
    ax.add_patch(
        patches.Rectangle(
            (1, min(atom_plddts)), xlim, 50-min(atom_plddts),
            color='#ff7c4f', alpha=1, label='Very Low (<50)'
        )
    )

    if len(chains := Counter(atom_chain_ids)) > 1:
        x = 1
        plt.axvline(x=x, color='grey', linestyle='--', alpha=1)
        ticks = []
        for chain_len in chains.values():
            ticks += [
                (x+chain_pos, str(chain_pos))
                for chain_pos in range(chain_len)
                if chain_pos % 1000 == 0
            ]
            plt.axvline(
                x=(x := x+chain_len), color='grey', linestyle='--', alpha=1
            )
        plt.xticks([t[0] for t in ticks], labels=[t[1] for t in ticks])
        plt.xlabel('Atom (per chain)', fontsize=12)
    else:
        plt.xlabel('Atom', fontsize=12)
    plt.ylabel('pLDDT', fontsize=12)
    plt.title(f'Predicted Local Distance Difference Test (pLDDT)', fontsize=14)
    plt.grid(alpha=0.5)

    ## as a scatter plot with each atom as a dot
    plt.plot(atom_plddts, color='black', marker='.', markersize=1, linewidth=0)
    ## as a line plot instead of scatter
    #plt.plot(atom_plddts, color='black', linewidth=1)

    plt.legend(
        title='Confidence Level', fontsize='small', bbox_to_anchor=(1.01, 1.0),
        loc='upper left'
    )
    plt.tight_layout()


def plot_pae(pae: list[int], token_chain_ids: list[str], show_chains= True) -> None:
    ''' Plot pae data '''

    plt.figure(figsize=(8, 8))
    plt.title(f'Predicted Aligned Error (PAE)', fontsize=14)
    plt.xlabel('Scored Residue', fontsize=12)
    plt.ylabel('Aligned Residue', fontsize=12)

    plt.imshow(np.array(pae), cmap='Greens_r')
    
    if show_chains:
        if len(chains := Counter(token_chain_ids)) > 1:
            xy = 0
            for v in list(chains.values())[:-1]:
                plt.axvline(x=(xy := xy+v), color='white', linestyle='--')
                plt.axhline(y=xy, color='white', linestyle='--')

    plt.colorbar(
        label='Expected Position Error (Ångströms)', orientation='horizontal',
        shrink=0.75, pad=0.075
    )
    plt.tight_layout()