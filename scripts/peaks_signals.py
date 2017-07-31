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

    print('Processing normalization by all library mapped reads')
    sizes = pd.read_csv(sizes_path, sep='\t', names=('name', 'sizes_pm'),
                        index_col='name')
    sizes_pm = sizes / 1000000
    coverage = pd.merge(coverage, sizes_pm, left_on="name", how='left',
                        right_index=True)

    print('Processing normalization by reads mapped to peaks')
    sizes_peaks_pm = coverage.groupby(["name"])\
        .agg({'coverage': 'sum'})\
        .rename(columns={'coverage': "sizes_peaks_pm"}) / 1000000
    coverage = pd.merge(coverage, sizes_peaks_pm, left_on="name", how='left',
                        right_index=True)

    print('Sizes RPM: {}'.format(sizes_pm))
    coverage['rpm'] = coverage['coverage'] / coverage['sizes_pm']
    save_signal(id, coverage, 'rpm', 'Saved RPM')

    print('Sizes peaks RPM: {}'.format(sizes_peaks_pm))
    coverage['rpm_peaks'] = coverage['coverage'] / coverage['sizes_peaks_pm']
    save_signal(id, coverage, 'rpm_peaks', 'Saved normalized reads by RPM reads in peaks signal')

    print('Processing RPKM normalization')
    coverage['rpk'] = (coverage['end'] - coverage['start']) / 1000.0
    coverage['rpkm'] = coverage['rpm'] / coverage['rpk']
    save_signal(id, coverage, 'rpkm', 'Saved RPKM')

    print('Processing RPKM_PEAKS normalization')
    coverage['rpkm_peaks'] = coverage['rpm_peaks'] / coverage['rpk']
    save_signal(id, coverage, 'rpkm_peaks', 'Saved normalized reads by RPKM reads in peaks signal')


def save_signal(id, coverage, signal_type, msg):
    pivot_df = pd.pivot_table(coverage, index=['chr', 'start', 'end'],
                              columns='name', values=signal_type, fill_value=0)
    result_filename = '{}_{}.tsv'.format(id, signal_type)
    pivot_df.to_csv(result_filename, sep='\t')
    print('{} to {}'.format(msg, result_filename))


def main():
    argv = sys.argv
    opts, args = getopt.getopt(argv[1:], "h", ["help"])
    if len(args) != 3:
        usage()
        sys.exit(1)

    process(args[0], args[1], args[2])


if __name__ == "__main__":
    main()
