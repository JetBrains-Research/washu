#!/usr/bin/env python
# Author: oleg.shpynov@jetbrains.com

import argparse
import os
from collections import namedtuple
import numpy as np
import pandas as pd
from bed.bedtrace import run, Bed

Record = namedtuple('Record', ['name', 'bdg'])


def process(regions, records, out):
    intersection_path = '{}.intersection.tsv'.format(out)
    if not os.path.exists(intersection_path):
        print('Save regions as BED3 format {}.bed3'.format(regions))
        regions3 = '{}.bed3'.format(regions)
        Bed(regions).save3(regions3)
        print('Compute summary intersection file {}'.format(intersection_path))
        cmd = [['bedtools', 'intersect', '-wa', '-wb', '-a', regions3, '-b', *[r.bdg for r in records], '-names',
                *[r.name for r in records], '-sorted'], ['awk', '-v', "OFS=\\t", '{print($1,$2,$3,$4,$8)}']]
        print(' | '.join([' '.join(c) for c in cmd]) + ' > ' + intersection_path)
        with open(intersection_path, "w") as intersection:
            run(cmd, stdout=intersection)
    print('Compute summary signal by {}'.format(intersection_path))
    compute_signal(intersection_path, out, records)


def compute_signal(intersection_path, out, records):
    intersection = pd.read_csv(intersection_path, names=('chr', 'start', 'end', 'name', 'coverage'), sep='\t')
    coverage = intersection.groupby(['chr', 'start', 'end', 'name'], as_index=False).aggregate(np.sum)
    # Pivot by names
    pivot = pd.pivot_table(coverage, index=['chr', 'start', 'end'], columns='name', values='coverage')
    raw_signal = '{}.csv'.format(out)
    pivot.to_csv(raw_signal)
    print('Saved raw signal to {}'.format(raw_signal))
    sizes = {}
    for r in records:
        print('Processing {}'.format(r))
        size = int(run([['cat', r.bdg], ['awk', '{cov+=$4} END{print(cov)}']])[0].decode("utf-8"))
        sizes[r] = size
        print('Size: {}'.format(size))
    print('TODO[shpynov] normalization by library coverage and by summary coverage of peaks')


def main():
    parser = argparse.ArgumentParser(description=
                                     'Process signal data of given bedgraph files within provided bed regions.')
    parser.add_argument('--regions', metavar='bed', type=str, help='BED file with regions of interest', required=True)
    parser.add_argument('--out', type=str, help='output prefix', required=True)
    parser.add_argument('--bedgraphs', metavar='name bdg', type=str, nargs='+',
                        help='Bedgraph files to be processed in format: name bdg [name bdg ...]', required=True)
    args = parser.parse_args()
    regions = args.regions
    print('Regions: {}'.format(regions))
    out = args.out
    print('Output suffix: {}'.format(out))
    # Load records
    it = iter(args.bedgraphs)
    records = [Record(*z) for z in zip(it, it)]
    print('Records: {}'.format(records))
    process(regions, records, out)


if __name__ == "__main__":
    main()
