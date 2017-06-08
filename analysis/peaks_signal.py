#!/usr/bin/env python
# Author: oleg.shpynov@jetbrains.com

import argparse
import os
from collections import namedtuple
import pandas as pd

from bed.bedtrace import run, Bed

Record = namedtuple('Record', ['name', 'bdg'])


def process(regions, records, out):
    sizes_path = '{}.sizes.csv'
    print('Processing sizes of libraries {}'.format(sizes_path))
    if not os.path.exists(sizes_path):
        with open(sizes_path, "w") as sizes_file:
            for r in records:
                size = int(run([['cat', r.bdg], ['awk', '{cov+=$4} END{print(cov)}']])[0].decode("utf-8"))
                print('{} size: {}'.format(r.bdg, size))
                sizes_file.write('{}\t{}\n'.format(r.name, size))

    intersection_path = '{}.intersection.tsv'.format(out)
    if not os.path.exists(intersection_path):
        regions3 = '{}.bed3'.format(regions)
        print('Save regions as BED3 format {}'.format(regions3))
        Bed(regions).save3(regions3)

        print('Compute summary intersection file {}'.format(intersection_path))
        cmd = [['bedtools', 'intersect', '-wa', '-wb', '-a', regions3, '-b', *[r.bdg for r in records], '-names',
                *[r.name for r in records], '-sorted'],
               ['awk', '-v', "OFS=\\t", '{print($1,$2,$3,$4,$8)}'],
               ['awk' "BEGIN{c="";s=0;e=0;n="";x=0} "
                "{ if ($1!=c || $2!=s || $3!=e || $4!=n) {"
                "if (x!=0) print($1,$2,$3,$4,x); c=$1;s=$2;e=$3;n=$4;x=$5 "
                "} else {x+=$5}} "
                "END {print($1,$2,$3,$4,x)}"]]
        print(' | '.join([' '.join(c) for c in cmd]) + ' > ' + intersection_path)
        with open(intersection_path, "w") as intersection_file:
            run(cmd, stdout=intersection_file)
    print('Compute summary signal by {}'.format(intersection_path))
    compute_signal(intersection_path, sizes_path, out)


def compute_signal(intersection_path, sizes_path, out):
    coverage = pd.read_csv(intersection_path, names=('chr', 'start', 'end', 'name', 'coverage'), sep='\t')
    # Pivot by names
    pivot = pd.pivot_table(coverage, index=['chr', 'start', 'end'], columns='name', values='coverage', fill_value=0)
    raw_signal = '{}.csv'.format(out)
    pivot.to_csv(raw_signal)
    print('Saved raw signal to {}'.format(raw_signal))

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
