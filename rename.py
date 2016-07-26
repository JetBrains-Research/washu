#!/usr/bin/env python
__author__ = 'oleg.shpynov@jetbrains.com'

import getopt
import sys
import os
import pandas as pd
import numpy as np
from subprocess import call

help_message = """
Script to rename samples to human readable names.
Usage: rename <folder> <go>?
"""


def usage():
    print(help_message)


def process(folder, go):
    """Process folder and renames files according to reads and samples tables."""
    if not go:
        print("Dry run, use `rename <folder> go` to rename.")
    reads = pd.read_csv('reads.csv')
    samples = pd.read_csv('samples.csv')
    for dirpath, dirs, files in os.walk(folder):
        print(dirpath)
        for f in files:
            n = f.replace('_s_', '_')
            records = reads[np.logical_and(
                reads['Run'].map(lambda x: x in n),
                np.logical_or(reads['ID'].map(lambda x: str(x) in n),
                              reads['TAG'].map(lambda x: len(str(x)) > 0 and str(x) in n)))]
            if len(records) != 1:
                print('Cannot rename', f)
            else:
                id = records.iloc[0]['ID']
                sample = samples.iloc[id - 1]['Sample']
                new_file = (str(id) + '_' + sample + '.fq')\
                    .replace(' ', '_').replace('(', '').replace(')', '').lower()
                print(f, new_file)
                if go:
                    call(['mv', dirpath + '/' + f, dirpath + '/' + new_file])


def main():
    argv = sys.argv
    opts, args = getopt.getopt(argv[1:], "h", ["help"])
    if len(args) == 0:
        usage()
        sys.exit(1)
    process(args[0], len(args) >= 2 and args[1] == 'go')


if __name__ == "__main__":
    main()
