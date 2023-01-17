#!/bin/bash
#title           :prep_database.sh
#description     :Prepare database for foldseek-fishing
#author          :Carlos Vigil Vásquez
#date            :20230117
#version         :20230117a
#notes           :Requires foldseek
#copyright       :Copyright (C) 2022 Carlos Vigil Vásquez (carlos.vigil.v@gmail.com)
#license         :Permission to copy and modify is granted under the MIT license

# Constants
GIT_ROOT=$(git rev-parse --show-toplevel)
DATABASE_PATH="${GIT_ROOT}/data/foldseek_dbs"

# ARgument parsing
if [[ -z $1 || $# -lt 2 || $# -ge 3 ]]; then
	cat <<-EOF
	Usage: ./${0##*/} database identifier
	
	Available databases in Foldseek:
	EOF
	foldseek databases | head -n9 | tail -n+3
	exit 1
fi

echo "prep_database.sh - $(date)"

# Download database
DATABASE="$1"
NAME="$2"

mkdir -p "$DATABASE_PATH"
cd "$DATABASE_PATH" || exit 1
foldseek databases "$DATABASE" "$NAME" /tmp

# Create FASTA files for aa and 3di sequences
foldseek convert2fasta "${NAME}" "${NAME}.fasta" 
foldseek lndb "${NAME}_h" "${NAME}_ss_h" 
foldseek convert2fasta "${NAME}_ss" "${NAME}_ss.fasta" 
