from pymol import cmd
import sys
import glob
from multiprocessing import Pool
import tqdm
import os


def pdb2cif(pdb):
    try:
        # Load structure
        cmd.load(pdb, os.path.basename(pdb))
        cmd.save(pdb.replace(".pdb", ".cif"), os.path.basename(pdb))
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
        for results in pool.map(pdb2cif, pdbs):
            pbar.update(1)
