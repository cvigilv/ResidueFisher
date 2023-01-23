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

    allpdbs = glob(source + '/*.pdb')
    for hit in sorted(allhits):
        print(hit)
        try:
            shutil.copy2(f"{source}/{hit[0:4]}.pdb", f"{target}/{hit[0:4]}.pdb")
        except FileNotFoundError:
            pattern = f"{source}/{hit[0:4].upper()}"
            closest = difflib.get_close_matches(pattern, allpdbs, n=1)[0]
            shutil.copy2(closest, f"{target}/{hit}.pdb")

if __name__ == "__main__":
    main()
