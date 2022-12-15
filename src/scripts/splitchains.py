from pymol import cmd
import sys
import glob
from multiprocessing import Pool
import tqdm
import os

pdbs = glob.glob(sys.argv[1])
print(len(pdbs))

def prepare_pdb(pdb):
    try:
        # Load structure
        cmd.load(pdb)

        # Remove useless objects
        cmd.remove("solvent")
        cmd.remove("hetatm")

        # Split chains
        cmd.split_chains()

        # Remove original structure
        if os.path.basename(pdb) == "cif":
            cmd.remove(pdb.strip(".cif"))
        else:
            cmd.remove(pdb.strip(".pdb"))

        # Save chains as new pdb structures
        for chain in list(cmd.get_object_list()):
            cmd.save("./" + chain + ".pdb", chain)
    except:
        print("ERROR:", pdb, "failed")
    finally:
        pass

    # Restart Pymol after processing
    cmd.reinitialize()

with Pool() as pool:
    pbar = tqdm.tqdm(total = len(pdbs))
    for results in pool.map(prepare_pdb, pdbs):
        pbar.update(1)
