#!/usr/bin/env python3

import sys
import numpy as np
from pymol import cmd
from Bio import SeqIO


def imprint(pdb_conservation_list, session_name):

    # Load and prepare structure
    for pdb, conservation in pdb_conservation_list:
        cmd.load(pdb)
        mol = pdb.strip(".pdb")
        obj = cmd.get_object_list(mol)[0]
        cmd.alter(mol, "b=0")  # Default conservation is zero

        # Add conservation information to structure
        bfacts = []
        residues = []
        # NOTE: This is to avoid problems with gaps
        for resi in [atom.resi for atom in cmd.get_model(("resi * and {}").format(mol)).atom]:
            if resi not in residues:
                residues.append(resi)

        for resn, line in enumerate(conservation):
            bfact = float(line)
            bfacts.append(bfact)
            cmd.alter(f"{mol} and resi {residues[resn]}", f"b={bfact}")

        # Construct conservation visualization
        cmd.show_as("cartoon", mol)
        cmd.cartoon("putty", mol)
        cmd.set("cartoon_putty_scale_min", 0.5, obj)

        cmd.color("tv_orange", "resn SER+THR+CYS+PRO+ASN+GLN")
        cmd.color("pink", "resn PHE+TYR+TRP")
        cmd.color("tv_green", "resn LYS+ARG+HIS")
        cmd.color("tv_blue", "resn GLY+ALA+VAL+LEU+MET+ILE")
        cmd.color("tv_red", "resn ASP+GLU")

        # # Highlight some interesting cases
        # if "_q" in mol:
        #     cmd.select(f"{mol} & byresidue b < 0")
        #     cmd.set_name("sele", f"conserved_gap_{mol}")

        if "_q" in mol:
            conservation_slices = [(0, 25), (25, 50), (50, 75), (75, 100)]
            for lower, upper in conservation_slices:
                cmd.select(f"{mol} & byres b > {lower} & byres b < {upper+0.1} ")
                cmd.set_name("sele", f"conservation_gt{lower}_le{upper}_{mol}")
            
            t_mol = "_".join(mol.replace("_q", "_t").split("_")[1:])
            cmd.select(f"byres *{t_mol} within 1.5 of conservation_gt75_le100_{mol}")
            cmd.set_name("sele", f"res_near_gt{lower}_le{upper}_{mol}")
    cmd.show("lines", f"conservation_gt75_le100*")
    cmd.deselect()
    cmd.save(session_name)

def id_sort(id):
    return int(id[0].split("_")[1]), int(id[0].split("_")[2]), id[0].split("_")[3]

def main():
    query = sys.argv[1]
    msa = sys.argv[2]
    consensus = sys.argv[3]
    conservation = sys.argv[4]
    output = sys.argv[5]

    query_aln = [(record.id, record.seq) for record in SeqIO.parse(msa, "fasta")]
    consensus_aln = [str(record.seq) for record in SeqIO.parse(consensus, "fasta")][0]
    pdb_conservation_list = []

    for frag_seqID in query_aln:
        if "fragment" in frag_seqID[0]:
            frag_pdb = ("{}.pdb").format(frag_seqID[0].split(":")[0])
            frag_seq = frag_seqID[1]
            with open(conservation) as io:
                io.readline()
                conservation_aln = (np.array(io.readline().split(",")).astype(np.float) * 100).astype(np.int)

            assert(len(frag_seq) == len(consensus_aln) == len(conservation_aln))

            cons_values = []
            for res in list(zip(frag_seq, consensus_aln, conservation_aln)):
                score = False
                if res[0] != '-':
                    score = res[2] if res[1] != '-' else -1

                if score:
                    cons_values.append(str(score))

            pdb_conservation = frag_pdb, cons_values
            pdb_conservation_list.append(pdb_conservation)
    pdb_conservation_list = sorted(pdb_conservation_list, key=id_sort, reverse=True)
    imprint(pdb_conservation_list, output)

main()
