#!/usr/bin/env python
# See scripts/peaks_signals.sh for data preprocessing
#
# Author: oleg.shpynov@jetbrains.com

import getopt
import sys

import numpy as np
import pandas as pd

help_message = 'ARGUMENTS: <coverage.csv> <sizes.csv> <id>'


def usage():
    print(help_message)


def process(coverage_path, sizes_path, id):
    print('PROCESSING ID', id)
    print('COVERAGE PATH', coverage_path)
    print('SIZES PATH', sizes_path)
    coverage = pd.read_csv(coverage_path, names=('chr', 'start', 'end', 'coverage', 'name'))

    print('Processing raw signal')
    pivot = pd.pivot_table(coverage, index=['chr', 'start', 'end'], columns='name', values='coverage', fill_value=0)
    raw_signal = '{}_raw.csv'.format(id)
    pivot.to_csv(raw_signal)
    print('Saved raw signal to {}'.format(raw_signal))

    print('Processing normalization by mln reads in library')
    sizes = pd.read_csv(sizes_path, names=('name', 'size'))
    sizes_mln = {row['name']: row['size'] / 1000000.0 for _, row in sizes.iterrows()}
    print('Sizes mln: {}'.format(sizes_mln))
    coverage['rpm'] = [row['coverage'] / sizes_mln[row['name']] for _, row in coverage.iterrows()]
    pivot_by_mln = pd.pivot_table(coverage, index=['chr', 'start', 'end'],
                                  columns='name', values='rpm', fill_value=0)
    raw_signal_by_mln = '{}_rpm.csv'.format(id)
    pivot_by_mln.to_csv(raw_signal_by_mln)
    print('Saved normalized reads per million RPM signal to {}'.format(raw_signal_by_mln))

    print('Processing normalization by mln reads in peaks')
    sizes_peaks_mln = {row['name']: np.sum(coverage[coverage['name'] == row['name']]['coverage']) / 1000000.0
                       for _, row in sizes.iterrows()}
    print('Sizes peaks mln: {}'.format(sizes_peaks_mln))
    coverage['rpm_peaks'] = [row['coverage'] / sizes_peaks_mln[row['name']] for _, row in
                                         coverage.iterrows()]
    pivot_by_coverage = pd.pivot_table(coverage, index=['chr', 'start', 'end'],
                                       columns='name', values='rpm_peaks', fill_value=0)
    raw_signal_by_coverage = '{}_rpm_peaks.csv'.format(id)
    pivot_by_coverage.to_csv(raw_signal_by_coverage)
    print('Saved normalized reads by mln reads in peaks signal to {}'.format(raw_signal_by_coverage))

    print('Processing RPKM normalization')
    coverage['coverage_rpkm'] = [row['coverage'] / ((row['end'] - row['start']) / 1000.0) / sizes_mln[row['name']] for
                                 _, row in coverage.iterrows()]
    pivot_rpkm = pd.pivot_table(coverage, index=['chr', 'start', 'end'],
                                columns='name', values='coverage_rpkm', fill_value=0)
    rpkm = '{}_rpkm.csv'.format(id)
    pivot_rpkm.to_csv(rpkm)
    print('Saved RPKM to {}'.format(rpkm))


def main():
    argv = sys.argv
    opts, args = getopt.getopt(argv[1:], "h", ["help"])
    if len(args) != 3:
        usage()
        sys.exit(1)

    process(args[0], args[1], args[2])


if __name__ == "__main__":
    main()
