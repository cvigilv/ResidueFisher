#!/usr/bin/env python3

import sys
import glob
import matplotlib.pyplot as plt

from pymol import cmd, stored
from collections import defaultdict

if __name__ == "__main__":
    print("fishingsummary.py")
    cases = glob.glob(sys.argv[1])
    kind = sys.argv[2]

    hit_counts = defaultdict(int)
    aa_counts = defaultdict(int)

    for session in cases:
        # Load session and identify kind of sequence used to "fish" residues
        hit = session.split("/")[-2].split("_")[1]

        # Get fished residues list
        cmd.load(session)
        cmd.select("resn in 'fished_*' and name CA")
        stored.fished_residues = []
        cmd.iterate("sele", "stored.fished_residues.append((resi, resn))")

        # Store fished residues per hit and fished residue frequency
        if stored.fished_residues != []:
            for resi, resn in stored.fished_residues:
                hit_counts[hit] += 1
                aa_counts[f"{resn}{resi}"] += 1

        # Close PyMOL session
        cmd.reinitialize()

    # Sort dictionaries
    hit_counts = dict(sorted(hit_counts.items(), key=lambda kv: kv[1], reverse=True))
    aa_counts = dict(sorted(aa_counts.items(), key=lambda kv: kv[1], reverse=True))

    # 'Fished' residues per hit {{{
    plt.figure(figsize=(0.25 * len(hit_counts), 4))

    plt.bar(list(hit_counts.keys()), list(hit_counts.values()), color="k")

    plt.xlabel("Hit")
    plt.ylabel("Fished residues")
    plt.xticks(rotation=90, ha="center")
    plt.title("'Fished' residues per hit")
    plt.legend(
        title="Sequence kind",
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        fancybox=True,
    )
    plt.savefig(
        sys.argv[3] + f"fished_residues.fished_per_hit.{kind}.png",
        dpi=300,
        bbox_inches="tight",
    ) # }}}
    # Frequency a residue is 'fished' {{{
    plt.figure(figsize=(0.25 * len(aa_counts), 4))
    plt.bar(list(aa_counts.keys()), list(aa_counts.values()), color="k")

    plt.xlabel("Query residue")
    plt.ylabel("Times fished")
    plt.xticks(rotation=90, ha="center")
    plt.title("Frequency a residue is 'fished'")
    plt.legend(
        title="Sequence kind",
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        fancybox=True,
    )
    plt.savefig(
        sys.argv[3] + f"fished_residues.fished_per_resi.{kind}.png",
        dpi=300,
        bbox_inches="tight",
    ) # }}}
