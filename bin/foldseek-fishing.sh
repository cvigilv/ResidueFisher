#!/bin/bash
#title           :foldseek-fishing.sh
#description     :Find structural and sequence information for PDB
#author          :Carlos Vigil Vásquez
#date            :20221214
#version         :20221214a
#notes           :Requires foldseek, moma2 docker, tmux and pymol
#copyright       :Copyright (C) 2022 Carlos Vigil Vásquez (cvigil2@uc.cl)
#license         :Permission to copy and modify is granted under the MIT license

# Constants
GIT_ROOT=$(git rev-parse --show-toplevel)
PDBS_PATH="${GIT_ROOT}/data/pdbs"
DATABASE_PATH="${GIT_ROOT}/data/foldseek_dbs"

if [[ -z $1 || $# -lt 2 || $# -ge 3 ]]; then
	cat <<-EOF
	Usage: ./${0##*/} pdb database
	EOF
	exit 1
fi

# 0. Initialize script
echo "foldseek-fishing.sh - $(date)"
QUERY=$(basename "$1" | sed "s/.pdb//" | sed "s/_//")
QUERY_PATH=$(readlink -f "$1")
DATABASE="$2"
RESULTS_PATH="${QUERY}_${DATABASE}"

mkdir "${GIT_ROOT}/results/${RESULTS_PATH}" || exit 1

# 1. Foldseek search
echo "1. Foldseek structural search"
FOLDSEEK_RESULTS="${GIT_ROOT}/results/${RESULTS_PATH}/foldseek"
mkdir "$FOLDSEEK_RESULTS" || exit 1
cd "$DATABASE_PATH" || exit 1

foldseek createdb "$QUERY_PATH" "$QUERY" --chain-name-mode 1 &> "$FOLDSEEK_RESULTS/alignment.log"
foldseek search "$QUERY" "$DATABASE" "$FOLDSEEK_RESULTS/alignment" /tmp \
  -e 0.001 \
  -s 6 \
  -k 6 \
  --tmscore-threshold 0 \
  -a \
  --max-seqs 100000 &> "$FOLDSEEK_RESULTS/alignment.log"
  # --cov-mode 2 \
  # -c 0.5 \
  # --max-seqs 100000 &> "$FOLDSEEK_RESULTS/alignment.log"
foldseek convertalis "$QUERY" "$DATABASE" \
  "$FOLDSEEK_RESULTS/alignment" \
  "$FOLDSEEK_RESULTS/alignment.tsv" \
  --format-output query,target,pident,alnlen,evalue,bits,qlen,tlen,qstart,qend,tstart,tend,qaln,taln,qseq,tseq,alntmscore &> "$FOLDSEEK_RESULTS/alignment.log"
cd - || exit 1


# 2. Foldseek analysis
echo "2. Foldseek structural search analysis"
mkdir "$FOLDSEEK_RESULTS/analysis" || exit 1
cd "$FOLDSEEK_RESULTS/analysis" || exit 1
python3 "${GIT_ROOT}/src/analyzefoldseek.py" "$FOLDSEEK_RESULTS/alignment.tsv" > "$FOLDSEEK_RESULTS/analysis/analysis.log"
cd - || exit 1

# 3. MSA
echo "3. MSA"
MSA_RESULTS="${GIT_ROOT}/results/${RESULTS_PATH}/msa"
mkdir "$MSA_RESULTS" || exit 1

# 3.1. MSA of aminoacid sequence
echo "3.1. Aminoacid-based MSA construction from Foldseek hits"
python3 "${GIT_ROOT}/src/getfasta.py" "$FOLDSEEK_RESULTS/alignment.tsv" "$DATABASE_PATH/${DATABASE}.fasta" "$MSA_RESULTS/hits.aa.fasta"
mafft "$MSA_RESULTS/hits.aa.fasta" 1> "$MSA_RESULTS/msa.aa.fasta" 2> "$MSA_RESULTS/msa.aa.log"
python3 "${GIT_ROOT}/src/analyzemsa.py" "$MSA_RESULTS/msa.aa.fasta" "$MSA_RESULTS/msa.aa_consensus.fasta"

# 3.2. MSA of structural sequence
echo "3.2. 3di-based MSA construction from Foldseek hits"
python3 "${GIT_ROOT}/src/getfasta.py" "$FOLDSEEK_RESULTS/alignment.tsv" "$DATABASE_PATH/${DATABASE}_ss.fasta" "$MSA_RESULTS/hits.3di.fasta"
mafft --aamatrix "$GIT_ROOT/src/3di.mat" "$MSA_RESULTS/hits.3di.fasta" 1> "$MSA_RESULTS/msa.3di.fasta" 2> "$MSA_RESULTS/msa.3di.log"
python3 "${GIT_ROOT}/src/analyzemsa.py" "$MSA_RESULTS/msa.3di.fasta" "$MSA_RESULTS/msa.3di_consensus.fasta"

# 3.3. Convert 3di to aminoacids, and viceversa
echo "3.3. Translation of 3di to aminoacid and viceversa"
python3 "${GIT_ROOT}/src/translatefasta.py" "$MSA_RESULTS/msa.3di.fasta" "$DATABASE_PATH/${DATABASE}.fasta"  "$MSA_RESULTS/msa.3di_aa.fasta"
python3 "${GIT_ROOT}/src/translatefasta.py" "$MSA_RESULTS/msa.aa.fasta" "$DATABASE_PATH/${DATABASE}_ss.fasta"  "$MSA_RESULTS/msa.aa_3di.fasta"
python3 "${GIT_ROOT}/src/analyzemsa.py" "$MSA_RESULTS/msa.3di_aa.fasta" "$MSA_RESULTS/msa.3di2aa_consensus.fasta"
python3 "${GIT_ROOT}/src/analyzemsa.py" "$MSA_RESULTS/msa.aa_3di.fasta" "$MSA_RESULTS/msa.aa23di_consensus.fasta"

# 4. Tree construction and filtering
echo "4. Tree construction and distance-based filtering"
TREE_RESULTS="${GIT_ROOT}/results/${RESULTS_PATH}/tree"
mkdir "$TREE_RESULTS" || exit 1

for SOURCE in aa 3di aa_3di 3di_aa;
do
	python3 "${GIT_ROOT}/src/createtree.py" "$MSA_RESULTS/msa.${SOURCE}.fasta" "$TREE_RESULTS/hits.${SOURCE}.nw" > "$TREE_RESULTS/tree.${SOURCE}.log"
	python3 "${GIT_ROOT}/src/filtertree.py" "$TREE_RESULTS/hits.${SOURCE}.nw" "$TREE_RESULTS/pruned_tree.${SOURCE}.nw" &> "$TREE_RESULTS/prunning.${SOURCE}.log"
done

# 5. Structural alignment
echo "5. Structural alignment of representatives"
mkdir -p ${GIT_ROOT}/results/${RESULTS_PATH}/moma/{input,output} || exit 1

echo "5.1. Download representative structures"
python3 "${GIT_ROOT}/src/getpdbs.py" "$TREE_RESULTS/pruned_tree.*" "$PDBS_PATH"

echo "5.2. MOMA alignment"
python3 "${GIT_ROOT}/src/movepdbs.py" "$TREE_RESULTS/pruned_tree.*" "$PDBS_PATH" "${GIT_ROOT}/results/${RESULTS_PATH}/moma/input"
tmux new-session -d -s "$QUERY"
tmux rename-window -t "$QUERY:1" 'moma'
tmux send-keys -t "$QUERY:1" "docker run -it -v $RESULTS_PATH:/home/momatools/data/ fggutierrez2018/moma2" C-m

# TODO: Implementar alineamiento estructural con MOMA2 (via tmux)
# TODO: Implementar imprinting de informacion de MSA a sesiones de Pymol

# Cleanup
rm ${DATABASE_PATH}/${QUERY}*
