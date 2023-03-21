#!/usr/bin/env python3

import sys
import glob
from pymol import cmd

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
    allhits = glob.glob(sys.argv[2] + "/*_t.pdb")
    output = sys.argv[3]

    # Load structures
    cmd.load(query)
    cmd.delete("conservation")
    query_name = cmd.get_names("objects")[0]

    # Stylize
    cmd.bg_color("white")
    cmd.set("antialias", 5)
    cmd.set("ambient", 0.7)
    cmd.set("ray_trace_mode", 1)

    # Fish residues
    for hit_file in allhits:
        # Get paths
        query_file = sys.argv[2] + "/" + query_name.replace("_","") + "_" + "_".join(hit_file.split("/")[-1].split("_")[1:-1]) + "_q.pdb"
        query_fragment = query_name.replace("_","") + "_" + "_".join(hit_file.split("/")[-1].split("_")[1:-1]) + "_q"
        hit_fragment = hit_file.strip('.pdb').split('/')[-1]
        fragment_name = "_".join(hit_file.split("/")[-1].split("_")[1:-1])

        # Load structures
        cmd.load(query_file)
        cmd.load(hit_file)

        # Generate fragment selections
        query_residues= f"({query_name} and conservation_gt75_le100 and name CA)"
        cmd.select(f"byres {query_fragment} within 0.1 of {query_residues}")
        cmd.set_name("sele", "conservation_gt75_le100_" + query_fragment)

        # Colors physicochemical groups
        for group, resn in AA_GROUPS.items():
            cmd.select(f"resn {resn} and {query_fragment} and name C*")
            cmd.set_name("sele", f"{group}_{query_fragment}")
            cmd.color(AA_COLORS[group], f"{group}_{query_fragment}")

            cmd.select(f"resn {resn} and {hit_fragment} and name C*")
            cmd.set_name("sele", f"{group}_{hit_fragment}")
            cmd.color(AA_COLORS[group], f"{group}_{hit_fragment}")

        # "Fish" residues
        for group in AA_GROUPS.keys():
            query_residues= f"({query_fragment} and conservation_gt75_le100_{query_fragment} and {group}_{query_fragment} and name CA)"
            hit_residues = f"({hit_fragment} and {group}_{hit_fragment} and name CA)"

            cmd.select(f"byres {query_residues} within {str(TOLERANCE)} of {hit_residues}")
            cmd.set_name("sele", f"fished_{group}_{fragment_name}")

            # Stylize fished residues
            cmd.show("lines", f"byres {hit_residues} within {str(TOLERANCE)} of fished_{group}_{fragment_name}")
            cmd.show("sticks", f"fished_{group}_{fragment_name}")

        # Remove junk selections
        cmd.delete(f"*_aa_{hit_fragment}")
        cmd.delete(f"*_aa_{query_fragment}")

        # Create scene for fragment for ease of analysis
        for fragment in [query_fragment, hit_fragment]:
            cmd.cartoon("loop", fragment)
            cmd.set("cartoon_transparency", 0.8, fragment)

        cmd.center(query_fragment)
        cmd.zoom(query_fragment)
        cmd.scene(f"Fragments_{fragment_name}", "store", f"Fished residues (conservation > 75% and same physicochemical properties) are shown in sticks.\nLines correspond to hit structure equivalent residues.")
        cmd.disable("all")

    # Create complete structure scenes
    cmd.enable("*_*_q")
    cmd.center("all")
    cmd.zoom("all")
    cmd.scene(f"Query", "store", f"Fished residues (conservation > 75% and same physicochemical properties) are shown in sticks.\nLines correspond to hit structure equivalent residues.")
    cmd.disable("all")

    cmd.enable("*_*_t")
    cmd.center("all")
    cmd.zoom("all")
    cmd.scene(f"Hit", "store", f"Fished residues (conservation > 75% and same physicochemical properties) are shown in sticks.\nLines correspond to hit structure equivalent residues.")
    cmd.disable("all")

    cmd.enable("*_*_*")
    cmd.center("all")
    cmd.zoom("all")
    cmd.scene(f"Complete", "store", f"Fished residues (conservation > 75% and same physicochemical properties) are shown in sticks.\nLines correspond to hit structure equivalent residues.")

    # Remove useless selections
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
