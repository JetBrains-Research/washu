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


def process(work_dir, id, sizes_path, peaks_sizes_path):
    data_path = os.path.join(work_dir, id, '{}.tsv'.format(id))
    data = pd.read_csv(data_path, sep='\t',
                       names=('chr', 'start', 'end', 'coverage', 'mean0', 'mean', 'name'))

    processing_chipseq = [n for n in data['name'] if re.match('.*input.*', n)] and \
                         not [n for n in data['name'] if re.match('.*(meth|trans|mirna).*', n)]

    if processing_chipseq:
        print("Processing ChIP-Seq, input found")
    else:
        print("Processing Methylation|Transcription|miRNA, input not found")

    print('Processing raw signal')
    pivot = pd.pivot_table(data, index=['chr', 'start', 'end'],
                           columns='name',
                           values='coverage' if processing_chipseq else 'mean',
                           fill_value=0)
    raw_signal = '{}_raw.tsv'.format(id)
    pivot.to_csv(raw_signal, sep='\t')
    print('Saved raw signal to {}'.format(raw_signal))

    # Normalization is available only for ChIP-Seq
    if not processing_chipseq:
        return

    print('Processing normalization by all library mapped reads')
    sizes = pd.read_csv(sizes_path, sep='\t', names=('name', 'sizes_pm'),
                        index_col='name')
    sizes_pm = sizes / 1000000
    data = pd.merge(data, sizes_pm, left_on="name", how='left', right_index=True)

    print('Sizes RPM: {}'.format(sizes_pm.to_string()))
    data['rpm'] = data['coverage'] / data['sizes_pm']
    save_signal(id, data, 'rpm', 'Saved RPM')

    print('Processing RPKM normalization')
    data['rpk'] = (data['end'] - data['start']) / 1000.0
    data['rpkm'] = data['rpm'] / data['rpk']
    save_signal(id, data, 'rpkm', 'Saved RPKM')

    if not(peaks_sizes_path and os.path.exists(peaks_sizes_path)):
        return
    print('Processing normalization by reads mapped to peaks')
    peaks_sizes = pd.read_csv(peaks_sizes_path, sep='\t', names=('name', 'sizes_peaks_pm'),
                              index_col='name')
    sizes_peaks_pm = peaks_sizes / 1000000
    print('Sizes peaks RPM: {}'.format(sizes_peaks_pm.to_string()))

    data = pd.merge(data, sizes_peaks_pm, left_on="name", how='left', right_index=True)
    data['rpm_peaks'] = data['coverage'] / data['sizes_peaks_pm']
    save_signal(id, data, 'rpm_peaks',
                'Saved normalized reads by RPM reads in peaks signal')
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
    if len(args) < 4:
        print("ARGUMENTS:  <work_dir> <id> <sizes.tsv> [<peaks.sizes.tsv>]\n"
              "     <work_dir>: folder with BW files\n"
              "     <sizes.tsv>: libraries sizes, used for RPM normalization\n"
              "     <peaks.sizes.tsv>: libraries sizes in peaks regions, used for RPM_peaks normalization\n"
              "CONVENTION: signal results are saved under <work_dir>/<id>\n\n"
              "ARGS: " + ",".join(args))
        sys.exit(1)

    work_dir = args[0]
    id = args[1]
    sizes_path = args[2]
    peaks_sizes_path = args[3] if len(args) == 4 else None
    print('WORK_DIR', work_dir)
    print('ID', id)
    print('SIZES PATH', sizes_path)
    print('PEAKS SIZES PATH', peaks_sizes_path)

    process(work_dir, id, sizes_path, peaks_sizes_path)

    signals_visualize.process(work_dir, id)


if __name__ == "__main__":
    main()
