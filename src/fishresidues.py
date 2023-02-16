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
    cmd.show_as("cartoon")
    cmd.cartoon("automatic")
    cmd.set("cartoon_transparency", 0.5, "all")
    for target in targets:
        query_file = sys.argv[2] + "/" + query_name.replace("_","") + "_" + "_".join(target.split("/")[-1].split("_")[1:-1]) + "_q.pdb"
        query_frag = query_name.replace("_","") + "_" + "_".join(target.split("/")[-1].split("_")[1:-1]) + "_q"
        target_frag = target.strip('.pdb').split('/')[-1]
        fragment_name = "_".join(target.split("/")[-1].split("_")[1:-1])
        cmd.load(query_file)
        cmd.load(target)
        query_residues1= f"({query_name} and conservation_gt50_le75 and name CA)"
        cmd.select(f"byres {query_frag} within 0.1 of {query_residues1}")
        cmd.set_name("sele", "conservation_gt50_le75_" + query_frag)
        query_residues2= f"({query_name} and conservation_gt75_le100 and name CA)"
        cmd.select(f"byres {query_frag} within 0.1 of {query_residues2}")
        cmd.set_name("sele", "conservation_gt75_le100_" + query_frag)
        # Colors and physico-chemical groups
        for group, resn in AA_GROUPS.items():
            cmd.select(f"resn {resn} and {query_frag}")
            cmd.set_name("sele", f"{group}_{query_frag}")
            cmd.color(AA_COLORS[group], f"{group}_{query_frag}")
            cmd.select(f"resn {resn} and {target_frag}")
            cmd.set_name("sele", f"{group}_{target_frag}")
            cmd.color(AA_COLORS[group], f"{group}_{target_frag}")
        # "Fish" residues
        for group in AA_GROUPS.keys():
            query_residues= f"({query_frag} and conservation_gt75_le100_{query_frag} and {group}_{query_frag} and name CA)"
            target_residues = f"({target_frag} and {group}_{target_frag} and name CA)"
            cmd.select(f"byres {target_residues} within {str(TOLERANCE)} of {query_residues}")
            cmd.set_name("sele", f"fished_{group}_{fragment_name}")
            cmd.show("lines", f"byres {query_residues} within {str(TOLERANCE)} of fished_{group}_{fragment_name}")
            cmd.show("sticks", f"fished_{group}_{fragment_name}")
            cmd.set("stick_transparency", 0.5, "all")
        cmd.delete(f"*_aa_{target_frag}")
        cmd.delete(f"*_aa_{query_frag}")
        cmd.scene(f"Fragments_{fragment_name}", "store", f"Residues that show more than 75% of conservation while maintaining physicochemical properties\nare shown in lines (query) and sticks (target).")
        cmd.disable("all")

    cmd.scene(f"Fragments_{fragment_name}", "recall")
    cmd.delete(query_name)
    cmd.delete("conservation_gt75_le100")
    cmd.delete("conservation_gt50_le75")
    cmd.delete("conservation_gt25_le50")
    cmd.delete("conservation_gt0_le25")
    cmd.delete("conserved_gap")

    # Save session
    cmd.deselect()
    cmd.save(output)

if __name__ == "__main__":
    main()
