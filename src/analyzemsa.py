#!/usr/bin/env python3

import os
import sys
import time
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO

AA_PALETTE = [
    "#C8C8C8",
    "#145AFF",
    "#00DCDC",
    "#E60A0A",
    "#E6E600",
    "#00DCDC",
    "#E60A0A",
    "#EBEBEB",
    "#8282D2",
    "#0F820F",
    "#0F820F",
    "#145AFF",
    "#E6E600",
    "#3232AA",
    "#DC9682",
    "#FA9600",
    "#FA9600",
    "#B45AB4",
    "#3232AA",
    "#0F820F",
    "#FFFFFF",
]
AA_LIST = "ARNDCQEGHILKMFPSTWYV-"


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
    conservation = getconservation(freqs, N)

    # Save consensus sequence
    consensus_seq = "".join(
        list(map(lambda x: dict(zip(range(21), AA_LIST))[x], consensus))
    )
    with open(sys.argv[2], "w+") as io:
        io.write(f"> Consensus\n{consensus_seq}")

    # Plot
    _, ax = plt.subplots(nrows=2, figsize=(20, 4))
    ax[0].bar(
        range(0, L), conservation, align="edge", linewidth=0, color="black", width=1.0
    )
    ax[0].set_ylabel("Conservation")
    ax[0].set_xlim(0, len(consensus))
    ax[0].set_xticks(
        np.ceil(np.linspace(0, len(consensus), 10)).astype(int)
    )
    ax[0].set_xticklabels(
        np.ceil(np.linspace(0, len(consensus), 10)).astype(int)
    )

    sns.heatmap(
        np.expand_dims(consensus, axis=0),
        cmap=AA_PALETTE,
        cbar_kws=dict(
            orientation="horizontal",
            ticks=np.linspace(0, 20, 22) + 0.5,
            label="Aminoacid",
        ),
        ax=ax[1],
    )
    ax[1].set_yticklabels(["Consensus"])
    ax[1].set_xticks(
        np.ceil(np.linspace(0, len(consensus), 10)).astype(int),
    )
    ax[1].set_xticklabels(
        np.ceil(np.linspace(0, len(consensus), 10)).astype(int),
    )
    ax[1].collections[0].colorbar.ax.set_xticklabels(AA_LIST + " ")
    ax[0].set_title(os.path.basename(sys.argv[1]))

    plt.savefig(f"{output_path}/msa.conservation_{str(time.time()).split('.')[0]}.pdf", dpi=300)


if __name__ == "__main__":
    main()
