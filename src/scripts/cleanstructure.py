from pymol import cmd
import sys
import glob
from multiprocessing import Pool
import tqdm


def cleanpdb(pdb):
    try:
        # Load structure
        cmd.load(pdb)

        # Remove useless objects
        cmd.remove("solvent")
        cmd.remove("hetatm")

        # Save chains as new pdb structures
        for chain in list(cmd.get_object_list()):
            cmd.save("./" + chain + "_clean.pdb", chain)
    except:
        print("ERROR:", pdb, "failed")
    finally:
        pass

    # Restart Pymol after processing
    cmd.reinitialize()

if __name__ == "__main__":
    pdbs = glob.glob(sys.argv[1])

    with Pool() as pool:
        pbar = tqdm.tqdm(total = len(pdbs))
        for results in pool.map(cleanpdb, pdbs):
            pbar.update(1)
