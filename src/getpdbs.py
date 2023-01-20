#!/usr/bin/env python3

import os
import sys
from glob import glob
from multiprocessing import Pool

import requests
from ete3 import Tree
from pymol import cmd
from requests.models import HTTPError

AF_DB_VERSION = "v4"


def cif2pdb(pdb):
    try:
        # Load structure
        cmd.load(pdb, os.path.basename(pdb))
        cmd.save(pdb.replace(".cif", ".pdb"), os.path.basename(pdb))
    except:
        print("ERROR:", pdb, "failed")
    finally:
        pass

    # Restart Pymol after processing
    cmd.reinitialize()


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

    # Restart PyMOL after processing
    cmd.reinitialize()


def downloadstructure(hit):
    availablepdbs = [
        os.path.basename(pdb).lower() for pdb in glob(sys.argv[2] + "/*.pdb")
    ]
    if hit[0:4].lower() + ".pdb" not in availablepdbs:
        iscif = False

        print(f"Downloading and processing {hit}")
        if "AF-" in hit:
            # Download PDB file of AlphaFold model
            hit = hit.strip(".gz")
            model_url = f"https://alphafold.ebi.ac.uk/files/{hit.strip('.gz')}"
            response = requests.get(model_url)
        else:
            # Download PDB or CIF file for hit structure
            try:
                hit = hit.split("_")[0] + ".pdb"
                response = requests.get(f"https://files.rcsb.org/download/{hit}")
                response.raise_for_status()
            except HTTPError:
                hit = hit.replace(".pdb", ".cif")
                response = requests.get(f"https://files.rcsb.org/download/{hit}")
                iscif = True

        with open(f"{sys.argv[2]}/{hit}", "wb") as f:
            f.write(response.content)

        # Convert CIF files to PDB spec
        if iscif:
            cif2pdb(f"{sys.argv[2]}/{hit}")

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
