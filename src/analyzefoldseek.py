#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
# import matplotlib as mpl
import numpy as np

import sys
# import os

def getstats(df, metrics = ['evalue', "alnlen", "pident", "alntmscore"]):
    return (
        df[metrics]
        .agg(["count", np.mean, np.std, np.min, np.median, np.max])
        .transpose()
    )

def main():
    print("analyzefoldseek.py")
    # Load results database and preprocess
    aln = pd.read_csv(
        sys.argv[1],
        sep='\t',
        names = "query,target,pident,alnlen,evalue,bits,qlen,tlen,qstart,qend,tstart,tend,qaln,taln,qseq,tseq,alntmscore".split(','),
    )
    aln["query"] = aln["query"].str.strip(".pdb").str.upper()
    aln["target"] = aln["target"].str.strip(".pdb").str.upper()

    print()
    print("Table. Hits statistic for relevant metrics")
    print(getstats(aln).to_markdown(floatfmt=".3g"))
    print()

    # Plot results
    aln["p-evalue"] = -np.log10(aln["evalue"])
    fig, axdict = plt.subplot_mosaic("ABC", figsize = (12,4))
    sns.scatterplot(
        x = "alnlen",
        y = "pident",
        hue = "p-evalue",
        data = aln,
        ax = axdict["A"],
        palette = "Spectral",
        legend = None
    )
    sns.scatterplot(
        x = "alnlen",
        y = "alntmscore",
        hue = "p-evalue",
        data = aln,
        ax = axdict["B"],
        palette = "Spectral",
    )
    sns.scatterplot(
        x = "pident",
        y = "alntmscore",
        hue = "p-evalue",
        data = aln,
        ax = axdict["C"],
        palette = "Spectral",
        legend = None
    )
    plt.tight_layout()
    plt.savefig("regressions.png", dpi = 300)

main()
