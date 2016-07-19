#!/usr/bin/env python
__author__ = 'oleg.shpynov@jetbrains.com'

import getopt
import sys
import os
import pandas as pd
import numpy as np
import re
from subprocess import call

help_message = 'Script to rename samples to human readable names.'


def usage():
    print(help_message)


def process(folder):
    """Process folder and renames files according to reads and samples tables."""
    reads = pd.read_csv('reads.csv')
    samples = pd.read_csv('samples.csv')
    for dirpath, dirs, files in os.walk(folder):
        print(dirpath)
        for f in files:
            n = f.replace('_s_', '_')
            records = reads[np.logical_and(
                reads['Run'].map(lambda x: x in n),
                np.logical_or(reads['ID'].map(lambda x: x in n),
                              reads['TAG'].map(lambda x: len(x) > 0 and x in n)))]
            if len(records) != 1:
                print('Cannot rename', f)
            else:
                r = records.iloc[0]
                sample = samples.iloc[r['ID'] - 1]['Sample']
                new_f = re.sub('^.*' + r['TAG'], sample, f).replace(' ', '_').lower()
                print(f, new_f)
                call(['mv', dirpath + '/' + f, dirpath + '/' + new_f])


def main():
    argv = sys.argv
    opts, args = getopt.getopt(argv[1:], "h", ["help"])
    if len(args) != 1:
        usage()
        sys.exit(1)
    process(args[0])


if __name__ == "__main__":
    main()
