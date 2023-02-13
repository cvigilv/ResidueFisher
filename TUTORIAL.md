# Tutorial

This tutorial intendes to showcase how to (1) prepare a query protein to study and the
target database to search for homology, (2) run the protocol and (3) interpret the results.

For this tutorial we will use que query protein PDB code 3F3P, nucleoporin pair Nup85-Seh1,
and the complete PDB dataset available directly from `foldseek`.

## Step 1. Preparation
### Query protein preparation
Unlike Foldseek, this protocol is intended to study a single protein chain; therefore, in
order to use `foldseek-fishing`, one must first extract this from its original PDB file.

In the `src/scripts` folder, there is a script called `splitchains.py`, which extracts all
the strings from a particular PDB file and saves them as separate files for use in
foldseek-fishing.

*Note*: the recommended way of preparing and storing all the structure files is creating a
new folder in data called `queries` and run the chain splitting script inside this folder.

Using our examle query structure, do the following:

```bash
# Ensure we have the conda environment activated
conda activate foldseek-fishing

mkdir data/queries
cd data/queries
wget https://files.rcsb.org/download/3F3P.pdb
python ../../src/scripts/splitchains.py 3F3P.pdb
```

From this example, a total of 13 should be found inside the `data/queries` folder: 1 for
the original structure (`3F3P.pdb`) and 12 corresponding to the chains A through L of 3F3P.

~~~bash
$ tree .
.
├── 3F3P.pdb
├── 3F3P_A.pdb
├── 3F3P_B.pdb
├── 3F3P_C.pdb
├── 3F3P_D.pdb
├── 3F3P_E.pdb
├── 3F3P_F.pdb
├── 3F3P_G.pdb
├── 3F3P_H.pdb
├── 3F3P_I.pdb
├── 3F3P_J.pdb
├── 3F3P_K.pdb
└── 3F3P_L.pdb
~~~

### Database preparation

To prepare a database using Foldseek, run the script `bin/prep_database.sh` script as
follows:

```bash
sh bin/prep_database.sh <FOLDSEEK-DATABASE-NAME> <INTERNAL-DATABASE-NAME>

# Using our example target database, do the following:
sh bin/prep_database.sh PDB mypdb
```
To see the available datasets, run `bin/prep_database.sh` without arguments.

*NOTE:* To prepare a database from PDB files, please refer to [foldseek tutorial](https://github.com/steineggerlab/foldseek#databases).
In order for `foldseek-fishing` to work correctly, user created databases must be inside a
directory named `data/foldseek_dbs` in the proyect root and must contain FASTA files for the
aminoacid sequence and 3di sequence, which can be created as follows:
```sh
foldseek convert2fasta <USER-DB-NAME> <USER-DB-NAME>.fasta
foldseek lndb <USER-DB-NAME>_h <USER-DB-NAME>_ss_h
foldseek convert2fasta <USER-DB-NAME>_ss <USER-DB-NAME>_ss.fasta
```

## Step 2. Running the protocol

To use foldseek-fishing, run the script `bin/foldseek-fishing.sh` script as follows:
```sh
sh bin/foldseek-fishing.sh <PDB-FILE> <INTERNAL-DATABASE-NAME>

# Example using the previously prepared protein and dataset
sh bin/foldseek-fishing.sh data/queries/3F3P_C.pdb mypdb
```

This will generate a folder in `results` with the following structure:
```
results/3F3PC_pdb/
├── foldseek/
├── msa/
├── moma/
└── tree/
```

Inside each subdirectory, log files and result can be found in order to analyse and study the protein used in the protocol.

