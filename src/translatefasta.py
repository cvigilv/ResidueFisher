#!/usr/bin/env python3

import sys
from textwrap import wrap
from Bio import SeqIO


def translateseq(source, target):
    transformed_str = ""

    i = 0
    for char in source:
        if char != "-":
            transformed_str += target[i]
            i += 1
        elif char == "-":
            transformed_str += "-"
        else:
            print("Possible error...")

    return transformed_str


def main():
    print("translatefasta.py")

    # Load fasta files
    fasta_file = sys.argv[1]
    target_dictionary = sys.argv[2]

    fasta = {record.id: record.seq for record in SeqIO.parse(fasta_file, "fasta")}
    dictionary = {
        record.id: record.seq for record in SeqIO.parse(target_dictionary, "fasta")
    }

    # Translate sequences
    translated_seqs = []
    for key, seq in fasta.items():
        translated_seq = translateseq(seq, dictionary[key])
        translated_seq = "\n".join(
            wrap(translated_seq, width=65, break_on_hyphens=False)
        )
        translated_seqs.append(f">{key}\n{translated_seq}")

    # Save
    with open(sys.argv[3], "w+") as io:
        io.write("\n".join(translated_seqs))


if __name__ == "__main__":
    main()
