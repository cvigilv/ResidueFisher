#!/usr/bin/env python3

import sys
import glob
import matplotlib.pyplot as plt

from pymol import cmd, stored
from collections import defaultdict

if __name__ == "__main__":
    print("fishingsummary.py")
    cases = glob.glob(sys.argv[1])
    print(cases)

    data = dict()
    fished_counts = defaultdict(int)

    for session in cases:
        hit = session.split("/")[-2].split("_")[1]
        cmd.load(session)
        cmd.select("resn in 'fished_*' and name CA")
        stored.fished_residues = []
        cmd.iterate("sele", "stored.fished_residues.append((resi, resn))")
        data[hit] = len(stored.fished_residues)
        for resi, resn in stored.fished_residues:
            fished_counts[f"{resn}{resi}"] += 1

        cmd.delete("(all)")

    data = dict(sorted(data.items(), key = lambda kv: kv[1], reverse = True))

    # Fished residues per query-hit pair
    plt.figure(figsize = (6,4))
    plt.bar(list(data.keys())[0:20], list(data.values())[0:20], color='k')
    plt.xlabel("Hit")
    plt.ylabel("Fished residues")
    plt.xticks(rotation=90, ha="center")
    plt.title("Top 20 hits with most 'fished' residues")
    plt.savefig(sys.argv[2] + "fished_residues.number_per_pair.png", dpi = 300, bbox_inches = "tight")


    # Fished residues counts
    fished_counts = dict(sorted(fished_counts.items(), key = lambda kv: kv[1], reverse = True))
    plt.figure(figsize = (6,4))
    plt.bar(list(fished_counts.keys())[0:20], list(fished_counts.values())[0:20],color='k')
    plt.xlabel("Query residue")
    plt.ylabel("Times fished")
    plt.xticks(rotation=90, ha="center")
    plt.title("Top 20 most 'fished' residues")
    plt.savefig(sys.argv[2] + "fished_residues.total_frequency.png", dpi = 300, bbox_inches = "tight")

