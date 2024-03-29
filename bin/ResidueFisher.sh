#!/bin/bash
#title           :ResidueFisher
#description     :Find structural and sequence information for PDB
#author          :Carlos Vigil Vásquez
#date            :20230324
#version         :1.0
#notes           :Requires foldseek, moma2 docker, tmux and pymol
#copyright       :Refer to LICENSE
#license         :Permission to copy and modify is granted under the MIT license

set -o errexit

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
echo "ResidueFisher - $(date) v1.0"
QUERY=$(basename "$1" | sed "s/.pdb//" | sed "s/_//")
QUERY_PATH=$(readlink -f "$1")
DATABASE="$2"
RESULTS_PATH="${QUERY}_${DATABASE}"

mkdir "${GIT_ROOT}/results/${RESULTS_PATH}" || exit 1

# 0. Foldseek search
echo "1. Foldseek structural search"
FOLDSEEK_RESULTS="${GIT_ROOT}/results/${RESULTS_PATH}/foldseek"
mkdir "$FOLDSEEK_RESULTS" || exit 1
cd "$DATABASE_PATH" || exit 1

# 1.1 Prepare query structure database
foldseek createdb "$QUERY_PATH" "$QUERY" --chain-name-mode 1 &> "$FOLDSEEK_RESULTS/alignment.log"
foldseek convert2fasta "${QUERY}" "${QUERY}.fasta" &> "$FOLDSEEK_RESULTS/alignment.log"
foldseek lndb "${QUERY}_h" "${QUERY}_ss_h" &> "$FOLDSEEK_RESULTS/alignment.log"
foldseek convert2fasta "${QUERY}_ss" "${QUERY}_ss.fasta" &> "$FOLDSEEK_RESULTS/alignment.log"

# 1.2 Structural search with Foldseek
foldseek search "$QUERY" "$DATABASE" "$FOLDSEEK_RESULTS/alignment" /tmp \
  -e 0.001 \
  -s 6 \
  -k 6 \
  --tmscore-threshold 0 \
  -a \
  --max-seqs 100000 &> "$FOLDSEEK_RESULTS/alignment.log"
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
python3 "${GIT_ROOT}/src/getfasta.py" "$FOLDSEEK_RESULTS/alignment.tsv" "$DATABASE_PATH/${DATABASE}.fasta" "/tmp/hits.aa.fasta"
cat "${GIT_ROOT}/data/foldseek_dbs/${QUERY}.fasta" "/tmp/hits.aa.fasta" > "$MSA_RESULTS/hits.aa.fasta"
mafft "$MSA_RESULTS/hits.aa.fasta" 1> "$MSA_RESULTS/msa.aa.fasta" 2> "$MSA_RESULTS/msa.aa.log"
python3 "${GIT_ROOT}/src/analyzemsa.py" "$MSA_RESULTS/msa.aa.fasta" "$MSA_RESULTS/msa.aa_consensus.fasta"

# 3.2. MSA of structural sequence
echo "3.2. 3di-based MSA construction from Foldseek hits"
python3 "${GIT_ROOT}/src/getfasta.py" "$FOLDSEEK_RESULTS/alignment.tsv" "$DATABASE_PATH/${DATABASE}_ss.fasta" "/tmp/hits.3di.fasta"
cat "${GIT_ROOT}/data/foldseek_dbs/${QUERY}_ss.fasta" "/tmp/hits.3di.fasta" > "$MSA_RESULTS/hits.3di.fasta"
mafft --aamatrix "$GIT_ROOT/src/3di.mat" "$MSA_RESULTS/hits.3di.fasta" 1> "$MSA_RESULTS/msa.3di.fasta" 2> "$MSA_RESULTS/msa.3di.log"
python3 "${GIT_ROOT}/src/analyzemsa.py" "$MSA_RESULTS/msa.3di.fasta" "$MSA_RESULTS/msa.3di_consensus.fasta"

# 3.4. Imprint MSA conservation to query structure
cd "$(dirname "$QUERY_PATH")"
for SOURCE in aa 3di;
do
	python3 "${GIT_ROOT}/src/conservationimprinting.py" \
		"$(basename "$QUERY_PATH")" \
		"$MSA_RESULTS/msa.$SOURCE.fasta" \
		"$MSA_RESULTS/msa.${SOURCE}_consensus.fasta" \
		"$MSA_RESULTS/msa.${SOURCE}_consensus.fasta_conservation" \
		"$MSA_RESULTS/${QUERY}.${SOURCE}_conservation.pse"
done
cd -

# 4. Tree construction and filtering
echo "4. Tree construction and distance-based filtering"
TREE_RESULTS="${GIT_ROOT}/results/${RESULTS_PATH}/tree"
mkdir "$TREE_RESULTS" || exit 1

for SOURCE in aa 3di;
do
	python3 "${GIT_ROOT}/src/createtree.py" "$MSA_RESULTS/msa.${SOURCE}.fasta" "$TREE_RESULTS/hits.${SOURCE}.nw" > "$TREE_RESULTS/tree.${SOURCE}.log"
	python3 "${GIT_ROOT}/src/filtertree.py" "$TREE_RESULTS/hits.${SOURCE}.nw" "$TREE_RESULTS/pruned_tree.${SOURCE}.nw" &> "$TREE_RESULTS/prunning.${SOURCE}.log"
done

# 5. Structural alignment
echo "5. Structural alignment of representatives"
mkdir -p "${GIT_ROOT}"/results/"${RESULTS_PATH}"/moma/{input,output} || exit 1

echo "5.1. Download representative structures"
python3 "${GIT_ROOT}/src/getpdbs.py" "$TREE_RESULTS/pruned_tree.*" "$PDBS_PATH"
python3 "${GIT_ROOT}/src/getpdb.py" "$QUERY" "$GIT_ROOT/results/$RESULTS_PATH/moma/input/"

echo "5.2. Collecting all structures"
python3 "${GIT_ROOT}/src/movepdbs.py" "$TREE_RESULTS/pruned_tree.*" "$PDBS_PATH" "$GIT_ROOT/results/$RESULTS_PATH/moma/input"

tmux new-session -d -s "$QUERY"
tmux send-keys -t "$QUERY:0" "docker run -it -v $GIT_ROOT/results/$RESULTS_PATH/moma/:/home/momatools/data/ fggutierrez2018/moma2" C-m
tmux send-keys -t "$QUERY:0" "cd /home/momatools/src" C-m

while read -r HIT_CHAIN; do
	HIT_PDB=$(echo "$HIT_CHAIN" | cut -d'_' -f1)
	HIT_CHAIN=$(echo "$HIT_CHAIN" | cut -d'_' -f2)
	QUERY_PDB=$(basename "${QUERY_PATH//.pdb/}"| cut -d'_' -f1)
	QUERY_CHAIN=$(basename "${QUERY_PATH//.pdb/}" | cut -d'_' -f2)
	PAIR="${QUERY_PDB}${QUERY_CHAIN}_${HIT_PDB}$HIT_CHAIN"
	echo "$PAIR"

	tmux send-keys -t "$QUERY:0" "python /home/momatools/src/MOMA2_pw.py -q /home/momatools/data/input/${QUERY_PDB}.pdb -t /home/momatools/data/input/${HIT_PDB}.pdb --cq $QUERY_CHAIN --ct $HIT_CHAIN -s both" C-m
	tmux send-keys -t "$QUERY:0" "python /home/momatools/src/generate_p1m.py /home/momatools/data/output/pairwise_alignments/${PAIR}_both" C-m
done<"$GIT_ROOT"/results/"$RESULTS_PATH"/moma/input/hits.txt
tmux send-keys -t "$QUERY:0" "> /home/momatools/data/output/.flag && exit" C-m
until [ -f "$GIT_ROOT/results/$RESULTS_PATH/moma/output/.flag" ]; do sleep 5; done
tmux kill-session -t "$QUERY"

echo "5.3. Fish interesting residues"
DIR="$GIT_ROOT"/results/"$RESULTS_PATH"/moma/output/
if [ "$(ls -A "$DIR")" ]; then
    mkdir "$GIT_ROOT/results/$RESULTS_PATH/moma/results/" || exit 1
    cd "$GIT_ROOT/results/$RESULTS_PATH/moma/results/"

    for i in "$GIT_ROOT"/results/"$RESULTS_PATH"/moma/output/pairwise_alignments/*_*_both;
    do
        PAIR=${i##*/}
        mkdir "$PAIR"
        cd "$PAIR"
        echo "Analyzing MOMA results in -> ${PAIR}"
        python3 "${GIT_ROOT}/src/fishresidues.py" \
            "${GIT_ROOT}/results/${RESULTS_PATH}/msa/${QUERY}.aa_conservation.pse" \
            "${GIT_ROOT}/results/${RESULTS_PATH}/moma/output/pairwise_alignments/$PAIR/best_combination" \
            "${QUERY}.aa_conservation_moma.pse"

        python3 "${GIT_ROOT}/src/fishresidues.py" \
            "${GIT_ROOT}/results/${RESULTS_PATH}/msa/${QUERY}.3di_conservation.pse" \
            "${GIT_ROOT}/results/${RESULTS_PATH}/moma/output/pairwise_alignments/$PAIR/best_combination" \
            "${QUERY}.3di_conservation_moma.pse"
        cd -
    done

    python3 "${GIT_ROOT}/src/fishingsummary.py" \
        "$GIT_ROOT/results/$RESULTS_PATH/moma/results/*/*aa*.pse" \
        "aa" \
        "$GIT_ROOT/results/$RESULTS_PATH/moma/"

    python3 "${GIT_ROOT}/src/fishingsummary.py" \
        "$GIT_ROOT/results/$RESULTS_PATH/moma/results/*/*3di*.pse" \
        "3di" \
        "$GIT_ROOT/results/$RESULTS_PATH/moma/"
else
    echo "Finishing foldseek-fishing because no alignments were found in MOMA"
fi

# Cleanup
rm "${DATABASE_PATH}"/"${QUERY}"*
