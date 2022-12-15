from pymol import cmd
import sys
import glob
from multiprocessing import Pool
import tqdm
import os


def cif2pdb(pdb):
    try:
        # Load structure
        cmd.load(pdb, os.path.basename(pdb))
        cmd.save(pdb.replace(".cif", ".pdb"), os.path.basename(pdb))
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
        for results in pool.map(cif2pdb, pdbs):
            pbar.update(1)
