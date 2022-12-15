#!/usr/bin/env python3

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import AlignIO

import sys

def main():
    print("createtree.py")

    aln = AlignIO.read(open(sys.argv[1]), 'fasta')
    print(aln)

    constructor = DistanceTreeConstructor()
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(aln)
    njtree = constructor.nj(dm)

    Phylo.write(njtree, sys.argv[2], "newick")

main()

