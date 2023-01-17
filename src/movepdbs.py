#!/usr/bin/env python3

import sys
import shutil

from ete3 import Tree
from glob import glob
import difflib

def insensitive_glob(pattern):
    def either(c):
        return '[%s%s]' % (c.lower(), c.upper()) if c.isalpha() else c
    return glob(''.join(map(either, pattern)))

def main():
    print("movepdbs.py")

    # Get representative hits for each tree
    trees = glob(sys.argv[1])
    source = sys.argv[2]
    target = sys.argv[3]

    allhits = set()
    for treepath in trees:
        treestring = open(treepath, "r").read().replace("\n", "")
        T = Tree(treestring, format=1)
        treehits = {leaf.name for leaf in T.get_leaves()}
        allhits.update(treehits)

    print(allhits)

    allpdbs = glob(source + '/*.pdb')
    for hit in allhits:
        pattern = f"{source}/{hit[0:4]}_{hit[5:].upper()}.pdb"
        # pattern2 = insensitive_glob(pattern)
        print(pattern, "->", difflib.get_close_matches(pattern, allpdbs))
        # shutil.copy2(pattern2[0], f"{target}/{hit}.pdb")


if __name__ == "__main__":
    main()
