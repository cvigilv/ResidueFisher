#!/usr/bin/env python3

import sys
import pandas as pd

from textwrap import wrap
from Bio import SeqIO

FOLDSEEK_HEADER = [
    "query",
    "target",
    "pident",
    "alnlen",
    "evalue",
    "bits",
    "qlen",
    "tlen",
    "qstart",
    "qend",
    "tstart",
    "tend",
    "qaln",
    "taln",
    "qseq",
    "tseq",
    "alntmscore",
]


def getsequence(seqid, start, end, fasta):
    seq = str(fasta[seqid])[int(start) - 1 : int(end) - 1]
    seq = "\n".join(wrap(seq, width=65, break_on_hyphens=False))

    return f"> {seqid}\n{seq}"


def main():
    print("getfasta.py")

    # Load results database and preprocess
    aln = pd.read_csv(
        sys.argv[1],
        sep="\t",
        names=FOLDSEEK_HEADER,
    )
    aln["query"] = aln["query"].apply(
        lambda x: x.split(".pdb")[0].lower()[::-1].capitalize()[::-1]
    )

    # Load sequences
    seqs = {record.id: record.seq for record in SeqIO.parse(sys.argv[2], "fasta")}

    # Create FASTA sequences per target
    fasta = (
        aln[["target", "tstart", "tend"]]
        .apply(lambda r: getsequence(r["target"], r["tstart"], r["tend"], seqs), axis=1)
        .values
    )

    # Save FASTA sequences
    with open(sys.argv[3], "w+") as io:
        io.write("\n".join(fasta))


if __name__ == "__main__":
    main()
