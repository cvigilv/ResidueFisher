#!/usr/bin/env python3

import shutil
import sys
from glob import glob

from ete3 import Tree


def insensitive_glob(pattern):
    def either(c):
        return "[%s%s]" % (c.lower(), c.upper()) if c.isalpha() else c

    return glob("".join(map(either, pattern)))


def main():
    print("movepdbs.py")

    # Get representative hits for each tree
    trees = glob(sys.argv[1])
    source_path = sys.argv[2]
    target_path = sys.argv[3]

    allhits = set()
    for treepath in trees:
        treestring = open(treepath, "r").read().replace("\n", "")
        T = Tree(treestring, format=1)
        treehits = {leaf.name for leaf in T.get_leaves()}
        allhits.update(treehits)

    with open(f"{target_path}/hits.txt", "w") as io:
        io.write("\n".join(allhits))

    for hit in allhits:
        shutil.copy2(f"{source_path}/{hit[0:4]}.pdb", f"{target_path}/{hit[0:4]}.pdb")

if __name__ == "__main__":
    main()
