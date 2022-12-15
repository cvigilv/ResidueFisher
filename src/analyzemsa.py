#!/usr/bin/env python3

import os
import sys
import time
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO

AA_PALETTE = [
    # Polar
    "#FF8C26",  # N Asn
    "#FF8C26",  # C Cys
    "#FF8C26",  # Q Gln
    "#FF8C26",  # P Pro
    "#FF8C26",  # S Ser
    "#FF8C26",  # T Thr
    # Non-Polar
    "#4D4DFF",  # A Ala
    "#4D4DFF",  # G Gly
    "#4D4DFF",  # I Ile
    "#4D4DFF",  # L Leu
    "#4D4DFF",  # M Met
    "#4D4DFF",  # V Val
    # Positive
    "#33FF33",  # R Arg
    "#33FF33",  # H His
    "#33FF33",  # K Lys
    # Negative
    "#FF3737",  # D Asp
    "#FF3737",  # E Glu
    # Aromatic
    "#FF7BFF",  # F Phe
    "#FF7BFF",  # W Trp
    "#FF7BFF",  # Y Tyr
    "#FFFFFF",  # - Gap
]
AA_LIST = "NCQPSTAGILMVRHKDEFWY-"
THRESHOLD = 0.3


def getaafreq(msa):
    L = len(msa[0].seq)
    N = len(msa)

    freq = np.zeros([L, 21])

    for i in range(0, N):
        for j in range(0, L):
            j_aa = AA_LIST.find(msa[i].seq[j])
            freq[j, j_aa] = freq[j, j_aa] + 1

    return freq


def getconsensus(freq):
    L = len(freq)
    consensus = np.zeros(L)
    for i in range(0, L):
        consensus[i] = freq[i].argmax()

    return consensus

def gethighlyconserved(freq):
    L = len(freq)
    N = sum(freq[0])
    consensus = np.zeros(L)
    for i in range(0, L):
        c = freq[i].argmax()
        consensus[i] = c if freq[i][c]/N > THRESHOLD else 20

    return consensus


def getconservation(freq, N):
    return np.sqrt(np.sum((np.square(freq / N - 0.05)), axis=1))


def main():
    print("analyzemsa.py")

    # Load results database and preprocess
    fasta_file = sys.argv[1]
    output_path = os.path.dirname(sys.argv[2])
    msa = list(SeqIO.parse(fasta_file, "fasta"))
    L = len(msa[0].seq)
    N = len(msa)

    # Calculate aa alignment metrics
    freqs = getaafreq(msa)
    consensus = getconsensus(freqs)
    highlyconserved = gethighlyconserved(freqs)
    conservation = getconservation(freqs, N)

    # Save consensus sequence
    consensus_seq = "".join(
        list(map(lambda x: dict(zip(range(21), AA_LIST))[x], consensus))
    )
    with open(sys.argv[2], "w+") as io:
        io.write(f"> Consensus\n{consensus_seq}")

    # Plot
    _, ax = plt.subplots(nrows=3, figsize=(int(L/100)*2, 6))
    ax[0].bar(
        range(0, L), conservation, align="edge", linewidth=0, color="black", width=1.0
    )
    ax[0].set_ylabel("Conservation")
    ax[0].set_xlim(0, len(consensus))
    ax[0].set_xticks(np.ceil(np.linspace(0, len(consensus), 10)).astype(int))
    ax[0].set_xticklabels(np.ceil(np.linspace(0, len(consensus), 10)).astype(int))
    ax[0].axhline(THRESHOLD, color="red", ls=":")
    ax[0].set_title(os.path.basename(sys.argv[1]))

    sns.heatmap(
        np.expand_dims(consensus, axis=0),
        cmap=AA_PALETTE,
        cbar=False,
        ax=ax[1],
        vmin=0,
        vmax=20,
    )
    ax[1].set_yticklabels(["Consensus"])
    ax[1].set_xticks(
        np.ceil(np.linspace(0, len(consensus), 10)).astype(int),
    )
    ax[1].set_xticklabels(
        np.ceil(np.linspace(0, len(consensus), 10)).astype(int),
    )

    sns.heatmap(
        np.expand_dims(highlyconserved, axis=0),
        cmap=AA_PALETTE,
        cbar_kws=dict(
            ticks=np.linspace(0, 20, 22) + 0.5,
            label="Aminoacid",
            use_gridspec=False,
            location="bottom"
        ),
        vmin=0,
        vmax=20,
        ax=ax[2],
    )
    ax[2].set_yticklabels(["Conserved"])
    ax[2].set_xticks(
       np.ceil(np.linspace(0, len(consensus), 10)).astype(int),
    )
    ax[2].set_xticklabels(
       np.ceil(np.linspace(0, len(consensus), 10)).astype(int),
    )
    ax[2].collections[0].colorbar.ax.set_xticklabels(AA_LIST + " ")

    # plt.tight_layout()
    plt.savefig(
        f"{output_path}/msa.conservation_{str(time.time()).split('.')[0]}.pdf", dpi=300
    )


if __name__ == "__main__":
    main()
