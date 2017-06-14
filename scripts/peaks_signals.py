#!/usr/bin/env python
# See scripts/peaks_signals.sh for data preprocessing
#
# Author: oleg.shpynov@jetbrains.com

import getopt
import sys

import numpy as np
import pandas as pd

help_message = 'ARGUMENTS: <coverage.tsv> <sizes.tsv> <id>'


def usage():
    print(help_message)


def process(coverage_path, sizes_path, id):
    print('PROCESSING ID', id)
    print('COVERAGE PATH', coverage_path)
    print('SIZES PATH', sizes_path)
    coverage = pd.read_csv(coverage_path, sep='\t', names=('chr', 'start', 'end', 'coverage', 'name'))

    print('Processing raw signal')
    pivot = pd.pivot_table(coverage, index=['chr', 'start', 'end'], columns='name', values='coverage', fill_value=0)
    raw_signal = '{}_raw.tsv'.format(id)
    pivot.to_csv(raw_signal, sep='\t')
    print('Saved raw signal to {}'.format(raw_signal))

    print('Processing normalization by RPM reads in library')
    sizes = pd.read_csv(sizes_path, sep='\t', names=('name', 'size'))
    sizes_rpm = {row['name']: row['size'] / 1000000.0 for _, row in sizes.iterrows()}
    print('Sizes RPM: {}'.format(sizes_rpm))
    coverage['rpm'] = [row['coverage'] / sizes_rpm[row['name']] for _, row in coverage.iterrows()]
    pivot_by_rpm = pd.pivot_table(coverage, index=['chr', 'start', 'end'],
                                  columns='name', values='rpm', fill_value=0)
    raw_signal_by_rpm = '{}_rpm.tsv'.format(id)
    pivot_by_rpm.to_csv(raw_signal_by_rpm, sep='\t', )
    print('Saved normalized reads per million RPM signal to {}'.format(raw_signal_by_rpm))

    print('Processing normalization by RPM reads in peaks')
    sizes_peaks_rpm = {row['name']: np.sum(coverage[coverage['name'] == row['name']]['coverage']) / 1000000.0
                       for _, row in sizes.iterrows()}
    print('Sizes peaks RPM: {}'.format(sizes_peaks_rpm))
    coverage['rpm_peaks'] = [row['coverage'] / sizes_peaks_rpm[row['name']] for _, row in
                                         coverage.iterrows()]
    pivot_by_coverage = pd.pivot_table(coverage, index=['chr', 'start', 'end'],
                                       columns='name', values='rpm_peaks', fill_value=0)
    raw_signal_by_coverage = '{}_rpm_peaks.tsv'.format(id)
    pivot_by_coverage.to_csv(raw_signal_by_coverage, sep='\t')
    print('Saved normalized reads by RPM reads in peaks signal to {}'.format(raw_signal_by_coverage))

    print('Processing RPKM normalization')
    coverage['coverage_rpkm'] = [row['coverage'] / ((row['end'] - row['start']) / 1000.0) / sizes_rpm[row['name']] for
                                 _, row in coverage.iterrows()]
    pivot_rpkm = pd.pivot_table(coverage, index=['chr', 'start', 'end'],
                                columns='name', values='coverage_rpkm', fill_value=0)
    rpkm = '{}_rpkm.tsv'.format(id)
    pivot_rpkm.to_csv(rpkm, sep='\t')
    print('Saved RPKM to {}'.format(rpkm))

    print('Sizes peaks RPKM: {}'.format(sizes_peaks_rpm))
    coverage['rpkm_peaks'] = [row['coverage'] / ((row['end'] - row['start']) / 1000.0) / sizes_peaks_rpm[row['name']] for
                                 _, row in coverage.iterrows()]
    pivot_by_coverage = pd.pivot_table(coverage, index=['chr', 'start', 'end'],
                                       columns='name', values='rpkm_peaks', fill_value=0)
    rpkm_peaks = '{}_rpkm_peaks.tsv'.format(id)
    pivot_by_coverage.to_csv(rpkm_peaks, sep='\t')
    print('Saved normalized reads by RPKM reads in peaks signal to {}'.format(rpkm_peaks))

def main():
    argv = sys.argv
    opts, args = getopt.getopt(argv[1:], "h", ["help"])
    if len(args) != 3:
        usage()
        sys.exit(1)

    process(args[0], args[1], args[2])


if __name__ == "__main__":
    main()
