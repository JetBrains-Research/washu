#!/usr/bin/env python
# See signals.sh for data preprocessing
#
# Author: oleg.shpynov@jetbrains.com

import math
import tempfile
import os
import subprocess
import pandas as pd

from scripts import signals_visualize
from scripts.util import *


def process(work_dir, id, sizes_path, peaks_sizes_path):
    data_path = os.path.join(work_dir, id, '{}.tsv'.format(id))
    data = pd.read_csv(data_path, sep='\t',
                       names=('chr', 'start', 'end', 'coverage', 'mean0', 'mean', 'name'))

    processing_chipseq = not [n for n in data['name'] if re.match('.*(meth|trans|mirna).*', n)]

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
    sizes = pd.read_csv(sizes_path, sep='\t', names=('name', 'size'), index_col='name')
    sizes_pm = sizes / 1000000
    sizes_pm.columns = ['size_pm']
    data = pd.merge(data, sizes_pm, left_on="name", how='left', right_index=True)

    data['rpm'] = data['coverage'] / data['size_pm']
    save_signal(id, data, 'rpm', 'Saved RPM')

    print('Processing RPKM normalization')
    data['rpk'] = (data['end'] - data['start']) / 1000.0
    data['rpkm'] = data['rpm'] / data['rpk']
    save_signal(id, data, 'rpkm', 'Saved RPKM')

    raw_data = pd.pivot_table(data, index=['chr', 'start', 'end'],
                              columns='name', values='coverage', fill_value=0)
    process_diffbind_scores(data_path, raw_data, sizes)

    if not (peaks_sizes_path and os.path.exists(peaks_sizes_path)):
        return
    print('Processing normalization by reads mapped to peaks')
    peaks_sizes = pd.read_csv(peaks_sizes_path, sep='\t', names=('name', 'sizes_peaks_pm'),
                              index_col='name')
    sizes_peaks_pm = peaks_sizes / 1000000

    data = pd.merge(data, sizes_peaks_pm, left_on="name", how='left', right_index=True)
    data['rpm_peaks'] = data['coverage'] / data['sizes_peaks_pm']
    save_signal(id, data, 'rpm_peaks',
                'Saved normalized reads by RPM reads in peaks signal')
    print('Processing RPKM_PEAKS normalization')
    data['rpkm_peaks'] = data['rpm_peaks'] / data['rpk']
    save_signal(id, data, 'rpkm_peaks',
                'Saved normalized reads by RPKM reads in peaks signal')


def save_signal(id, data, signal_type, msg):
    pivot_df = pd.pivot_table(data, index=['chr', 'start', 'end'],
                              columns='name', values=signal_type, fill_value=0)
    result_filename = '{}_{}.tsv'.format(id, signal_type)
    pivot_df.to_csv(result_filename, sep='\t')
    print('{} to {}'.format(msg, result_filename))


TMM_R_PATH = os.path.dirname(os.path.realpath(__file__)) + '/../R/tmm.R'


def process_diffbind_scores(data_path, data, sizes):
    print('Processing DiffBind scores')
    od_input = [c for c in data.columns.values if is_od_input(c)][0]
    yd_input = [c for c in data.columns.values if is_yd_input(c)][0]
    ods = [c for c in data.columns if is_od(c)]
    yds = [c for c in data.columns if is_yd(c)]
    records = [(d, od_input, OLD) for d in ods] + [(d, yd_input, YOUNG) for d in yds]
    scores, lib_sizes = process_scores(data, sizes, records)
    save_scores(data, scores, re.sub('.tsv', '_scores.tsv', data_path), 'diffbind scores')

    print('TMM normalization')
    scores_tmpfile = tempfile.NamedTemporaryFile(prefix='scores', suffix='.tsv').name
    scores.to_csv(scores_tmpfile, index=False, sep='\t')
    print('Saved scores to', scores_tmpfile)
    tmm_file = scores_tmpfile.replace('.tsv', '_tmm.tsv')
    sizes_tmpfile = tempfile.NamedTemporaryFile(prefix='sizes', suffix='.tsv').name
    lib_sizes.to_csv(sizes_tmpfile, index=False, sep='\t', header=None)
    print('Saved sizes to', sizes_tmpfile)
    print('TMM normalization using R')
    cmd = "Rscript " + TMM_R_PATH + " " + scores_tmpfile + " " + sizes_tmpfile + " " + tmm_file
    subprocess.run(cmd, shell=True)
    # Difference between DBA_SCORE_TMM_MINUS_FULL and DBA_SCORE_TMM_MINUS_FULL_CPM is in bCMP
    print('TMM Scores DBA_SCORE_TMM_MINUS_FULL_CPM')
    scores_tmm = pd.read_csv(tmm_file, sep='\t') * 10000000
    save_scores(data, scores_tmm, re.sub('.tsv', '_scores_tmm.tsv', data_path), 'diffbind TMM scores')


def save_scores(data, scores, file, name):
    scores.index = data.index
    scores.to_csv(file, sep='\t')
    print('{} to {}'.format(name, file))


def score(cond, cont, scale):
    """
Computes DiffBind score
See documents on how to compute scores
https://docs.google.com/document/d/1zH5cw5Zal546xkoFFCVqhhYmf3742efhddz5cqpD9PQ/edit?usp=sharing
"""
    if scale > 1:
        scale = 1
    if scale != 0:
        cont = math.ceil(cont * scale)
    # According to pv.get_reads() function, if reads number is < 1 than it should be 1
    return max(1, cond - cont)


def process_scores(data, sizes, records):
    scores_processed = pd.DataFrame()
    sizes_processed = pd.DataFrame(columns=['name', 'size'])
    for cond, cont, g in records:
        scale = sizes.loc[cond]['size'] / sizes.loc[cont]['size']
        prefix = '' if g is None else g.prefix
        scores_processed[prefix + cond] = [score(z[0], z[1], scale) for z in zip(data[cond], data[cont])]
        sizes_processed.loc[len(sizes_processed)] = (prefix + cond, sizes.loc[cond]['size'])
    return scores_processed, sizes_processed


def main():
    argv = sys.argv
    opts, args = getopt.getopt(argv[1:], "h", ["help"])
    if len(args) < 3:
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
