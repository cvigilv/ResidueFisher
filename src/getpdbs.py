#!/usr/bin/env python3

import os
import sys
import requests

from ete3 import Tree
from multiprocessing import Pool
from glob import glob
from pymol import cmd
from requests.models import HTTPError

AF_DB_VERSION = "v4"


def splitchains(pdb):
    # Load structure
    cmd.load(pdb)

    # Remove useless objects
    cmd.remove("solvent")
    cmd.remove("HETATM")

    # Split chains
    org_pdb = cmd.get_object_list()[0]
    cmd.split_chains()

    # Remove original structure
    cmd.remove(org_pdb)

    # Save chains as new pdb structures
    for chain in list(cmd.get_object_list()):
        cmd.save(os.path.dirname(pdb) + "/" + chain + ".pdb", chain)

    # Restart Pymol after processing
    cmd.reinitialize()


def downloadstructure(hit):
    run_splitting = False
    availablepdbs = [os.path.basename(pdb).strip(".pdb").strip(".cif") for pdb in glob(sys.argv[2] + "/*")]
    print(availablepdbs, list(glob(sys.argv[2]+"/*")))
    if hit not in availablepdbs or f"{hit[0:5]+hit[5:].upper()}" not in availablepdbs:
        print(f"Downloading and processing {hit}")
        if "AF-" in hit:
            hit = hit.strip(".gz")
            model_url = f"https://alphafold.ebi.ac.uk/files/{hit.strip('.gz')}"
            response = requests.get(model_url)
        else:
            run_splitting = True

            try:
                hit = hit.split("_")[0] + ".pdb"
                model_url = f"https://files.rcsb.org/download/{hit}"
                response = requests.get(model_url)
                response.raise_for_status()
            except HTTPError:
                hit = hit.replace(".pdb", ".cif")
                model_url = f"https://files.rcsb.org/download/{hit}"
                response = requests.get(model_url)

        with open(f"{sys.argv[2]}/{hit}", "wb") as f:
            f.write(response.content)

        if run_splitting:
            splitchains(f"{sys.argv[2]}/{hit}")
    else:
        print(f"{hit} already downloaded and processed")


def main():
    print("getpdbs.py")

    # Get representative hits for each tree
    trees = glob(sys.argv[1])

    allhits = set()
    for treepath in trees:
        treestring = open(treepath, "r").read().replace("\n", "")
        T = Tree(treestring, format=1)
        treehits = {leaf.name for leaf in T.get_leaves()}
        allhits.update(treehits)

    print(f"Number of representative hits: {len(allhits)}")

    with Pool() as p:
        _ = p.map(downloadstructure, allhits)


if __name__ == "__main__":
    main()
