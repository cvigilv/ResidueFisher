#!/usr/bin/env python3

import sys

from getpdbs import downloadstructure

def main():
    print("getpdb.py")

    id = sys.argv[1]
    if "AF-" in id:
        downloadstructure(id)
    else:
        downloadstructure(id[0:4])

if __name__ == "__main__":
    main()
