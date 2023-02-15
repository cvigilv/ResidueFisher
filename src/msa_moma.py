#!/usr/bin/env python3

import os
import sys
from Bio import SeqIO
import fnmatch
import warnings

def add_fragments_msa(PDBFile):
    warnings.simplefilter('ignore')
    with open(PDBFile, 'r') as pdb_file:
        for record in SeqIO.parse(pdb_file, 'pdb-atom'):
            header = ("{}:{}_fragment").format(PDBFile.split("/")[-1].strip(".pdb"), record.id.split(":")[1])
            return f">{header}\n{(record.seq).replace('X', '')}\n"

def main():
    moma_alignments_path = sys.argv[1]
    msa_alignments_path = sys.argv[2]
    source = sys.argv[3]
    for folder_query in os.listdir(moma_alignments_path):
        if fnmatch.fnmatch(folder_query, '*both'):
            if not os.path.exists(folder_query):
                os.makedirs(folder_query)
            fragments_seq_list = []
            folder_query_moma = ("{}{}/best_combination").format(moma_alignments_path,folder_query)
            with open(f"{folder_query}/fragments_to_add_{source}.fasta", "w") as outfile:
                for query in os.listdir(folder_query_moma):
                    if fnmatch.fnmatch(query, '*_*.pdb'):
                        file_path = ("{}/{}").format(folder_query_moma, query)
                        os.system(("cp {} {}").format(file_path, folder_query))
                        fragment_seq = add_fragments_msa(file_path)
                        outfile.write(fragment_seq)

        cmd_mafft = ("mafft --auto --addfragments {f}/fragments_to_add_{s}.fasta --keeplength {path}msa.{s}.fasta 1> {f}/frag_aln_{s}.fasta 2> {f}/frag_aln_{s}.log").format(f=folder_query, path=msa_alignments_path, s=source)
        os.system(cmd_mafft)
        cmd_cp_files = ("cp {path}msa.{s}_consensus.fasta {path}msa.{s}_consensus.fasta_conservation {f}").format(f=folder_query, path=msa_alignments_path, s=source)
        os.system(cmd_cp_files)

if __name__ == "__main__":
    main()
