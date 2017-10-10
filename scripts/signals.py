#!/usr/bin/env python
# See signals.sh for data preprocessing
#
# Author: oleg.shpynov@jetbrains.com

import getopt
import sys

import os
import pandas as pd
import re

from scripts import signals_visualize

help_message = 'ARGUMENTS: <coverage.tsv> <sizes.tsv> <id>'


def usage():
    print(help_message)


def process(data_path, sizes_path, id):
    print('PROCESSING ID', id)
    print('COVERAGE PATH', data_path)
    print('SIZES PATH', sizes_path)
    data = pd.read_csv(data_path, sep='\t',
                       names=('chr', 'start', 'end', 'coverage', 'mean0', 'mean', 'name'))

    chipseq_processing = [n for n in data['name'] if re.match('.*input.*', n)] and \
                         not [n for n in data['name'] if re.match('.*(meth|trans|mirna).*', n)]

    print('Processing raw signal')
    pivot = pd.pivot_table(data, index=['chr', 'start', 'end'],
                           columns='name',
                           values='coverage' if chipseq_processing else 'mean',
                           fill_value=0)
    raw_signal = '{}_raw.tsv'.format(id)
    pivot.to_csv(raw_signal, sep='\t')
    print('Saved raw signal to {}'.format(raw_signal))

    if not chipseq_processing:
        print("Processing METH|TRANS|MIRNA")
        return

    print('Processing normalization by all library mapped reads')
    sizes = pd.read_csv(sizes_path, sep='\t', names=('name', 'sizes_pm'),
                        index_col='name')
    sizes_pm = sizes / 1000000
    data = pd.merge(data, sizes_pm, left_on="name", how='left',
                    right_index=True)

    print('Processing normalization by reads mapped to peaks')
    sizes_peaks_pm = data.groupby(["name"]).agg({'coverage': 'sum'}).rename(
        columns={'coverage': "sizes_peaks_pm"}) / 1000000
    data = pd.merge(data, sizes_peaks_pm, left_on="name", how='left',
                    right_index=True)

    print('Sizes RPM: {}'.format(sizes_pm))
    data['rpm'] = data['coverage'] / data['sizes_pm']
    save_signal(id, data, 'rpm', 'Saved RPM')

    print('Sizes peaks RPM: {}'.format(sizes_peaks_pm))
    data['rpm_peaks'] = data['coverage'] / data['sizes_peaks_pm']
    save_signal(id, data, 'rpm_peaks',
                'Saved normalized reads by RPM reads in peaks signal')

    print('Processing RPKM normalization')
    data['rpk'] = (data['end'] - data['start']) / 1000.0
    data['rpkm'] = data['rpm'] / data['rpk']
    save_signal(id, data, 'rpkm', 'Saved RPKM')

    print('Processing RPKM_PEAKS normalization')
    data['rpkm_peaks'] = data['rpm_peaks'] / data['rpk']
    save_signal(id, data, 'rpkm_peaks',
                'Saved normalized reads by RPKM reads in peaks signal')


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

    coverage_path = args[0]
    sizes_path = args[1]
    id = args[2]
    print('Processing bed_signal.py {} {} {}'.format(coverage_path, sizes_path, id))
    process(coverage_path, sizes_path, id)

    bam_bw_folder = os.path.dirname(os.path.dirname(coverage_path))
    signals_visualize.process(bam_bw_folder, id)


if __name__ == "__main__":
    main()
