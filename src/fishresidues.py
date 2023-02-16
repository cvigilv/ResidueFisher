#!/usr/bin/env python3

import sys
from types import TracebackType
import numpy as np
import glob
from pymol import cmd, stored
from Bio import SeqIO

# Aminoacid classification
AA_GROUPS = {
    "polar_aa": "SER+THR+CYS+PRO+ASN+GLN",
    "negative_aa": "ASP+GLU",
    "aromatic_aa": "PHE+TYR+TRP",
    "positive_aa": "LYS+ARG+HIS",
    "nonpolar_aa": "GLY+ALA+VAL+LEU+MET+ILE",
}
AA_COLORS = {
    "polar_aa": "0xFF8C26",
    "nonpolar_aa": "0x4D4DFF",
    "positive_aa": "0x33FF33",
    "negative_aa": "0xFF3737",
    "aromatic_aa": "0xFF7BFF",
}
TOLERANCE = 1.5

def main():
    query = sys.argv[1]
    targets = glob.glob(sys.argv[2] + "/*_t.pdb")
    output = sys.argv[3]

    # Load structures
    cmd.load(query)
    cmd.delete("conservation")
    query_name = cmd.get_names("objects")[0]
    for target in targets:
        cmd.load(target)

    # Generate selections and color
    for group, resn in AA_GROUPS.items():
        cmd.select("resn " + resn)
        cmd.set_name("sele", group)
        cmd.color(AA_COLORS[group], group)

    # "Fish" residues
    for group in AA_GROUPS.keys():
        query_residues= f"({query_name} and conservation_gt75_le100 and {group} and name CA)"
        target_residues = f"(*_t and {group} and name CA)"
        cmd.select(f"byres {query_residues} within {str(TOLERANCE)} of {target_residues}")
        cmd.set_name("sele", f"fished_{group}")

    cmd.select(" + ".join([f"fished_{group}" for group in AA_GROUPS.keys()]))
    cmd.set_name("sele", "fished_aa")

    # Set style
    cmd.show("lines", "conservation_gt75_le100")
    cmd.show("sticks", "fished_aa")
    cmd.set("cartoon_transparency", 0.5, "all")
    cmd.bg_color("white")

    # Save session
    cmd.deselect()
    cmd.save(output)

if __name__ == "__main__":
    main()
