#!/usr/bin/env python
# Author: oleg.shpynov@jetbrains.com

import math
import tempfile

import numpy as np
import pandas as pd

from downstream.aging import *
from downstream.signals import signals_tests
from downstream.signals import signals_visualize
from scripts.util import *


def process(data_path, sizes_path, peaks_sizes_path, post_process_callback):
    loaded = pd.read_csv(data_path, sep='\t',
                         names=('chr', 'start', 'end', 'coverage', 'mean0', 'mean', 'name'))

    processing_chipseq = not [n for n in loaded['name'] if re.match('.*(meth|trans|mirna).*', n)]
    if processing_chipseq:
        print("Processing ChIP-Seq")
    else:
        print("Processing Methylation|Transcription|miRNA, input not found")

    print('Processing RAW signal')
    raw_data = pd.pivot_table(loaded, index=['chr', 'start', 'end'],
                              columns='name',
                              values='coverage' if processing_chipseq else 'mean',
                              fill_value=0)
    raw_path = re.sub('.tsv', '_raw.tsv', data_path)
    raw_data.to_csv(raw_path, sep='\t')
    print('Saved RAW signal to {}'.format(raw_path))
    post_process_callback(raw_path)

    print("Processing QUANTILE normalization")
    quantile_path = re.sub('.tsv', '_q.tsv', data_path)
    process_quantile(quantile_path, raw_data)
    post_process_callback(quantile_path)

    # Normalization is available only for ChIP-Seq
    if not processing_chipseq:
        return

    print('Processing RPM normalization')
    sizes = pd.read_csv(sizes_path, sep='\t', names=('name', 'size'), index_col='name')
    sizes_pm = sizes / 1000000.0
    sizes_pm.columns = ['size_pm']
    data = pd.merge(loaded, sizes_pm, left_on="name", how='left', right_index=True)

    data['rpm'] = data['coverage'] / data['size_pm']
    rpm_path = re.sub('.tsv', '_rpm.tsv', data_path)
    save_signal(rpm_path, data, 'rpm', 'Saved RPM')
    post_process_callback(rpm_path)

    print('Processing RPKM normalization')
    data['rpk'] = (data['end'] - data['start']) / 1000.0
    data['rpkm'] = data['rpm'] / data['rpk']
    rpkm_path = re.sub('.tsv', '_rpkm.tsv', data_path)
    save_signal(rpkm_path, data, 'rpkm', 'Saved RPKM')
    post_process_callback(rpkm_path)

    if peaks_sizes_path and os.path.exists(peaks_sizes_path):
        print('Processing NORM')
        raw_data = pd.pivot_table(loaded, index=['chr', 'start', 'end'],
                                  columns='name', values='coverage', fill_value=0)
        norm_path = re.sub('.tsv', '_norm.tsv', data_path)
        process_norm(norm_path, raw_data, sizes_path, peaks_sizes_path)
        post_process_callback(norm_path)

    print('Processing DiffBind scores')
    raw_data = pd.pivot_table(loaded, index=['chr', 'start', 'end'],
                              columns='name', values='coverage', fill_value=0)
    process_diffbind_scores(data_path, raw_data, sizes, post_process_callback)


def save_signal(path, data, signal_type, msg):
    pivot_df = pd.pivot_table(data, index=['chr', 'start', 'end'],
                              columns='name', values=signal_type, fill_value=0)
    pivot_df.to_csv(path, sep='\t')
    print('{} to {}'.format(msg, path))


def quantile_normalize_using_target(x, target):
    """Both `x` and `target` are numpy arrays of equal lengths."""
    target_sorted = np.sort(target)
    return target_sorted[x.argsort().argsort()]


def process_quantile(output, data):
    # Normalize everything to the first track
    signal_columns = [c for c in data.columns if not is_input(c)]
    for c in signal_columns[1:]:
        data[c] = quantile_normalize_using_target(data[c],
                                                  data[signal_columns[0]])
    data.to_csv(output, sep='\t')
    print('{} to {}'.format('Saved QUANTILE', output))


def process_norm(output, data, sizes_path, peaks_sizes_path):
    """Normalization on library depths, Frips, and intersection fraction with peaks"""
    sizes = pd.read_csv(sizes_path, sep='\t', names=('name', 'tags'), index_col='name')
    peaks_sizes = pd.read_csv(peaks_sizes_path, sep='\t', names=('name', 'tags_in_peaks'),
                              index_col='name')

    counts = pd.merge(sizes, peaks_sizes, left_index=True, right_index=True)

    # Here we assume that we have single input for everything
    input_column = [c for c in data.columns if is_input(c)][0]
    counts['input_tags'] = counts.loc[input_column, 'tags']
    counts['input_tags_in_peaks'] = counts.loc[input_column, 'tags_in_peaks']

    counts = counts[~ counts.index.str.contains("input")]

    no_peaks_coverages = counts['tags'] - counts['tags_in_peaks']
    no_peaks_input_coverages = counts['input_tags'] - counts['input_tags_in_peaks']

    input_coef = no_peaks_coverages / no_peaks_input_coverages
    counts['input_coef'] = input_coef

    counts['data'] = counts['tags_in_peaks'] - input_coef * counts['input_tags_in_peaks']

    mean_count = np.mean(counts['data'])

    not_input_columns = [c for c in data.columns if not is_input(c)]
    for column in not_input_columns:
        v = data[column] - input_coef[column] * data[input_column]
        data[column] = np.maximum(0, v * mean_count / counts['data'][column])

    data[not_input_columns].to_csv(output, sep='\t')
    print('{} to {}'.format('Saved norm', output))


TMM_R_PATH = os.path.dirname(os.path.realpath(__file__)) + '/tmm.R'


def process_diffbind_scores(data_path, data, sizes, post_process_callback):
    od_input = [c for c in data.columns if is_od_input(c)][0]
    yd_input = [c for c in data.columns if is_yd_input(c)][0]
    ods = [c for c in data.columns if is_od(c)]
    yds = [c for c in data.columns if is_yd(c)]
    records = [(d, od_input, OLD) for d in ods] + [(d, yd_input, YOUNG) for d in yds]
    scores, lib_sizes = process_scores(data, sizes, records)
    scores_path = re.sub('.tsv', '_scores.tsv', data_path)
    save_scores(scores_path, data, scores, 'diffbind scores')
    post_process_callback(scores_path)

    tmm_results = process_tmm(scores, lib_sizes)
    scores_tmm = tmm_results * 10000000
    tmm_scores_path = re.sub('.tsv', '_scores_tmm.tsv', data_path)
    save_scores(tmm_scores_path, data, scores_tmm, 'diffbind TMM scores')
    post_process_callback(tmm_scores_path)


def process_tmm(data, lib_sizes):
    print('Scores TMM normalization')
    scores_tmpfile = tempfile.NamedTemporaryFile(prefix='scores', suffix='.tsv').name
    data.to_csv(scores_tmpfile, index=False, sep='\t')
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
    tmm_results = pd.read_csv(tmm_file, sep='\t')
    return tmm_results


def save_scores(path, data, scores, name):
    scores.index = data.index
    scores.to_csv(path, sep='\t')
    print('{} to {}'.format(name, path))


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
    # According to DiffBind pv.get_reads() function (utils.R:239)
    # if reads number is < 1 than it should be 1
    return max(1, cond - cont)


def process_scores(data, sizes, records):
    scores_processed = pd.DataFrame()
    sizes_processed = pd.DataFrame(columns=['name', 'size'])
    for cond, cont, g in records:
        scale = sizes.loc[cond]['size'] / sizes.loc[cont]['size']
        prefix = '' if g is None else g.prefix
        scores_processed[prefix + cond] = \
            [score(z[0], z[1], scale) for z in zip(data[cond], data[cont])]
        sizes_processed.loc[len(sizes_processed)] = (prefix + cond, sizes.loc[cond]['size'])
    return scores_processed, sizes_processed


def main():
    argv = sys.argv
    opts, args = getopt.getopt(argv[1:], "h", ["help"])
    if len(args) < 2:
        print("ARGUMENTS:  <data.tsv> <sizes.tsv> [<peaks.sizes.tsv>]\n"
              "     <data.tsv>: tsv file from signal_bw.sh\n"
              "     <sizes.tsv>: libraries sizes, used for RPM normalization\n"
              "     <peaks.sizes.tsv>: libraries sizes in peaks"
              "                        for peaks coverage normalization\n\n"
              "ARGS: " + ",".join(args))
        sys.exit(1)

    data_path = args[0]
    sizes_path = args[1]
    peaks_sizes_path = args[2] if len(args) == 3 else None
    print('DATA PATH', data_path)
    print('LIBRARIES SIZES PATH', sizes_path)
    print('PEAKS SIZES PATH', peaks_sizes_path)

    def post_process_callback(path):
        print('RESULT', path)
        signals_visualize.process(path)
        signals_tests.process(path)

    process(data_path, sizes_path, peaks_sizes_path, post_process_callback)


if __name__ == "__main__":
    main()
