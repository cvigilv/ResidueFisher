#!/usr/bin/env python3

import sys
import numpy as np
from pymol import cmd
from Bio import SeqIO


def imprint(pdb, conservation, session_name):
    # Constants
    PALETTE = ["deeppurple", "white", "tv_orange"]

    # Load and prepare structure
    cmd.load(pdb)
    mol = pdb.strip(".pdb")
    obj = cmd.get_object_list(mol)[0]
    cmd.alter(mol, "b=0")  # Default conservation is zero

    # Add conservation information to structure
    bfacts = []
    residues = [
        int(atom.resi) for atom in cmd.get_model("name CA").atom
    ]  # This is to avoid problems with gaps
    for resn, line in enumerate(conservation):
        bfact = float(line)
        bfacts.append(bfact)
        cmd.alter(f"{mol} and resi {residues[resn]}", f"b={bfact}")

    # Construct conservation visualization
    cmd.show_as("cartoon", mol)
    cmd.cartoon("putty", mol)
    cmd.set("cartoon_putty_scale_min", 0.5, obj)
    cmd.spectrum(
        "b",
        " ".join(PALETTE),
        "%s" % mol,
        minimum=0,
        maximum=100,
    )
    cmd.ramp_new("conservation", obj, [0, 50, 100], PALETTE)
    cmd.recolor()

    # Highlight some interesting cases
    cmd.select(f"{mol} & byresidue b < 0")
    cmd.set_name("sele", "conserved_gap")

    conservation_slices = [(0, 25), (25, 50), (50, 75), (75, 100)]
    for lower, upper in conservation_slices:
        cmd.select(f"{mol} & byresidue b > {lower} & byresidue b < {upper+0.1} ")
        cmd.set_name("sele", f"conservation_gt{lower}_le{upper}")

    cmd.show("lines", "conservation_gt75_le100")
    cmd.deselect()
    cmd.save(session_name)


def main():
    query = sys.argv[1]
    msa = sys.argv[2]
    consensus = sys.argv[3]
    conservation = sys.argv[4]
    output = sys.argv[5]

    query_aln = [(record.id, record.seq) for record in SeqIO.parse(msa, "fasta")][0][1]
    consensus_aln = [str(record.seq) for record in SeqIO.parse(consensus, "fasta")][0]
    with open(conservation) as io:
        io.readline()
        conservation_aln = (np.array(io.readline().split(",")).astype(np.float) * 100).astype(np.int)

    assert(len(query_aln) == len(consensus_aln) == len(conservation_aln))

    cons_values = []
    for res in list(zip(query_aln, consensus_aln, conservation_aln)):
        score = False
        if res[0] != '-':
            score = res[2] if res[1] != '-' else -1

        if score:
            cons_values.append(str(score))

    # print(cons_values)

    imprint(query, cons_values, output)

main()
