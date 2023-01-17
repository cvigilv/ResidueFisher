# Foldseek-fishing

`foldseek-fishing` is a bioinformatics protocol for the search of protein homology using a three-step ”search, detect, and enrich” model that uses a combination of structural and sequence aligners working in tandem to filter and enrich conservation signals.

## Table of Contents

- [Dependencies](#dependencies)
- [Installation](#installation)
- [Usage](#usage)
- [Support](#support)
- [Contributing](#contributing)

## Dependencies

`foldseek-fishing` depends on:
- conda
- docker
- mafft
- tmux

## Installation

To install `foldseek-fishing`, run the following code snippet:
```sh
git clone https://github.com/cvigilv/foldseek-fishing
cd foldseek-fishing
make configure
```
In order to use `foldseek-fishing`, que conda environment must be active (`conda activate foldseek-fishing`)

## Usage
### Database preparation

To prepare a database using Foldseek, run the script `bin/prep_database.sh` script as follows:

```sh
sh bin/prep_database.sh <FOLDSEEK-DATABASE-NAME> <INTERNAL-DATABASE-NAME>

# Example
sh bin/prep_database.sh PDB mypdb
```

To prepare a database from PDB files, please refer to [foldseek tutorial](https://github.com/steineggerlab/foldseek#databases).

In order for `foldseek-fishing` to work correctly, databases must be inside a directory named `data/foldseek_dbs` in the proyect root.

### Foldseek-fishing Usage

To use foldseek-fishing, run the script `bin/foldseek-fishing.sh` script as follows:
```sh
sh bin/foldseek-fishing.sh <PDB-FILE> <INTERNAL-DATABASE-NAME>

# Example
sh bin/foldseek-fishing.sh 3F3P_C.pdb mypdb
```

This will generate a folder in `results` with the following structure:
```
results/3F3PC_pdb/
├── foldseek/
├── msa/
└── tree/
```

Inside each subdirectory, log files and result can be found in order to analyse and study the protein used in the protocol.

## Support

Please [open an issue](https://github.com/cvigilv/foldseek-fishing/issues/new) for
support.

## Contributing

Please contribute using [Github Flow]
(https://guides.github.com/introduction/flow/). Create a branch, add
commits, and [open a pull request](https://github.com/cvigilv/foldseek-fishing/compare/).

## License

MIT

