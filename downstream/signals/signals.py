#!/usr/bin/env python
# Author: oleg.shpynov@jetbrains.com

import math
import multiprocessing
import tempfile

import numpy as np
import pandas as pd
from sklearn import preprocessing

from downstream.aging import *
from downstream.signals import signals_tests
from downstream.signals import signals_visualize
from scripts.util import *


def process(data_path, sizes_path, peaks_sizes_path, *, processed=4):
    pool = multiprocessing.Pool(processes=processed)
    loaded = pd.read_csv(data_path, sep='\t',
                         names=('chr', 'start', 'end', 'coverage', 'mean0', 'mean', 'name'))

    processing_chipseq = not [n for n in loaded['name'] if re.match('.*(meth|trans|mirna).*', n)]
    if processing_chipseq:
        print("Processing ChIP-Seq")
    else:
        print("Processing Methylation|Transcription|miRNA, input not found")
    pool.apply_async(raw_normalization,
                     args=(data_path, loaded, processing_chipseq, post_process_callback),
                     error_callback=error_callback)

    if processing_chipseq:
        pool.apply_async(rpm_normalization,
                         args=(data_path, loaded, sizes_path, post_process_callback),
                         error_callback=error_callback)

        if peaks_sizes_path and os.path.exists(peaks_sizes_path):
            pool.apply_async(frip_normalization,
                             args=(data_path, loaded, sizes_path,
                                   peaks_sizes_path, post_process_callback),
                             error_callback=error_callback)

        pool.apply_async(diffbind_normalization,
                         args=(data_path, loaded, sizes_path, post_process_callback),
                         error_callback=error_callback)

    pool.close()
    pool.join()


def raw_normalization(data_path, loaded, processing_chipseq, post_process_callback):
    """Raw signal with Quantile normalization and standard scaling"""
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
    quantile_path = re.sub('.tsv', '_rawq.tsv', data_path)
    process_quantile(quantile_path, raw_data)
    post_process_callback(quantile_path)

    print("Processing Z normalization")
    z_path = re.sub('.tsv', '_rawz.tsv', data_path)
    process_z(z_path, raw_data)
    post_process_callback(z_path)


def rpm_normalization(data_path, loaded, sizes_path, post_process_callback):
    """RPM/RPKM normalization"""
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


def save_signal(path, data, signal_type, msg):
    pivot_df = pd.pivot_table(data, index=['chr', 'start', 'end'],
                              columns='name', values=signal_type, fill_value=0)
    pivot_df.to_csv(path, sep='\t')
    print('{} to {}'.format(msg, path))


def process_z(output, data):
    """Center to the mean and component wise scale to unit variance."""
    z_df = pd.DataFrame()
    for c in [c for c in data.columns if not is_input(c)]:
        z_df[c] = preprocessing.scale(data[c])
    z_df.index = data.index
    z_df.to_csv(output, sep='\t')
    print('{} to {}'.format('Saved Z', output))


def quantile_normalize_using_target(x, target):
    """Both `x` and `target` are numpy arrays of equal lengths."""
    target_sorted = np.sort(target)
    return target_sorted[x.argsort().argsort()]


def process_quantile(output, data):
    quantile_df = pd.DataFrame()
    # Normalize everything to the first track
    signal_columns = [c for c in data.columns if not is_input(c)]
    target_column = signal_columns[0]
    target = data[target_column]
    print("Quantile normalization to", target_column)
    quantile_df[target_column] = target
    for c in signal_columns[1:]:
        quantile_df[c] = quantile_normalize_using_target(data[c], target)
    quantile_df.to_csv(output, sep='\t')
    print('{} to {}'.format('Saved QUANTILE', output))


def frip_normalization(data_path, loaded, sizes_path, peaks_sizes_path, post_process_callback):
    print('Processing FRIP normalization')
    data = pd.pivot_table(loaded, index=['chr', 'start', 'end'],
                          columns='name', values='coverage', fill_value=0)
    """Normalization on library depths, FRIPs"""
    sizes = pd.read_csv(sizes_path, sep='\t', names=('name', 'tags'), index_col='name')
    peaks_sizes = pd.read_csv(peaks_sizes_path, sep='\t', names=('name', 'tags_in_peaks'),
                              index_col='name')

    counts = pd.merge(sizes, peaks_sizes, left_index=True, right_index=True)

    # Here we assume that we have single input for everything
    input_column = [c for c in data.columns if is_input(c)][0]
    not_input_columns = [c for c in data.columns if not is_input(c)]
    counts['input_tags'] = counts.loc[input_column, 'tags']
    counts['input_tags_in_peaks'] = counts.loc[input_column, 'tags_in_peaks']

    # Omit input columns
    counts = counts[~ counts.index.str.contains("input")]

    no_peaks_coverages = counts['tags'] - counts['tags_in_peaks']
    no_peaks_input_coverages = counts['input_tags'] - counts['input_tags_in_peaks']

    # Input coefficient shows the ratio of signal track noise to input track level.
    # The smaller coefficient the better ChIP-Seq tracks we deal with.
    input_coef = no_peaks_coverages / no_peaks_input_coverages
    counts['input_coef'] = input_coef
    counts['data'] = counts['tags_in_peaks'] - input_coef * counts['input_tags_in_peaks']

    # Subtract number number of noisy reads
    frip_df = pd.DataFrame()
    for column in not_input_columns:
        frip_df[column] = np.maximum(0, data[column] - input_coef[column] * data[input_column])

    frip_path = re.sub('.tsv', '_frip.tsv', data_path)
    frip_df[not_input_columns].to_csv(frip_path, sep='\t')
    print('{} to {}'.format('Saved FRIP normalized', frip_path))
    post_process_callback(frip_path)

    # Per million reads normalization
    frip_pm = pd.DataFrame()
    for column in not_input_columns:
        frip_pm[column] = frip_df[column] * 1000000.0 / counts['data'][column]
    frip_linear_path = re.sub('.tsv', '_fripm.tsv', data_path)
    frip_pm[not_input_columns].to_csv(frip_linear_path, sep='\t')
    print('{} to {}'.format('Saved FRIP per million reads normalized', frip_linear_path))
    post_process_callback(frip_linear_path)

    # Quantile normalization
    frip_quantile_path = re.sub('.tsv', '_fripq.tsv', data_path)
    process_quantile(frip_quantile_path, data)
    post_process_callback(frip_quantile_path)

    # Z normalization
    frip_z_path = re.sub('.tsv', '_fripz.tsv', data_path)
    process_z(frip_z_path, data)
    post_process_callback(frip_z_path)


def diffbind_normalization(data_path, loaded, sizes_path, post_process_callback):
    print('Processing DiffBind scores')
    data = pd.pivot_table(loaded, index=['chr', 'start', 'end'],
                          columns='name', values='coverage', fill_value=0)
    sizes = pd.read_csv(sizes_path, sep='\t', names=('name', 'size'), index_col='name')

    od_input = [c for c in data.columns if is_od_input(c)][0]
    yd_input = [c for c in data.columns if is_yd_input(c)][0]
    ods = [c for c in data.columns if is_od(c)]
    yds = [c for c in data.columns if is_yd(c)]
    records = [(d, od_input, OLD) for d in ods] + [(d, yd_input, YOUNG) for d in yds]
    scores, lib_sizes = diffbind_scores(data, sizes, records)
    scores_path = re.sub('.tsv', '_scores.tsv', data_path)
    scores.index = data.index
    scores.to_csv(scores_path, sep='\t')
    print('Saved Diffbind scores to {}'.format(scores_path))
    post_process_callback(scores_path)

    # Quantile scores normalization
    scores_quantile_path = re.sub('.tsv', '_scoresq.tsv', data_path)
    process_quantile(scores_quantile_path, scores)
    post_process_callback(scores_quantile_path)

    # Z normalization
    scores_z_path = re.sub('.tsv', '_scoresz.tsv', data_path)
    process_z(scores_z_path, scores)
    post_process_callback(scores_z_path)

    # TMM normalization as in original DiffBind
    tmm_results = process_tmm(scores, lib_sizes)
    scores_tmm = tmm_results * 10000000
    scores_tmm_path = re.sub('.tsv', '_scores_tmm.tsv', data_path)
    scores_tmm.index = data.index
    scores_tmm.to_csv(scores_path, sep='\t')
    print('Saved Diffbind TMM scores to {}'.format(scores_tmm_path))
    post_process_callback(scores_tmm_path)


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
    cmd = "Rscript " + os.path.dirname(os.path.realpath(__file__)) + '/tmm.R' + \
          " " + scores_tmpfile + " " + sizes_tmpfile + " " + tmm_file
    subprocess.run(cmd, shell=True)
    # Difference between DBA_SCORE_TMM_MINUS_FULL and DBA_SCORE_TMM_MINUS_FULL_CPM is in bCMP
    print('TMM Scores DBA_SCORE_TMM_MINUS_FULL_CPM')
    tmm_results = pd.read_csv(tmm_file, sep='\t')
    return tmm_results


def diffbind_scores(data, sizes, records):
    """
    Computes DiffBind score
    See documents on how to compute scores
    https://docs.google.com/document/d/1zH5cw5Zal546xkoFFCVqhhYmf3742efhddz5cqpD9PQ/edit?usp=sharing
    """

    def score(condition, control, scale):
        if scale > 1:
            scale = 1
        if scale != 0:
            control = math.ceil(control * scale)
        # According to DiffBind pv.get_reads() function (utils.R:239)
        # if reads number is < 1 than it should be 1
        return max(1, condition - control)

    scores = pd.DataFrame()
    sizes_processed = pd.DataFrame(columns=['name', 'size'])
    for cond, cont, g in records:
        s = sizes.loc[cond]['size'] / sizes.loc[cont]['size']
        prefix = '' if g is None else g.prefix
        scores[prefix + cond] = \
            [score(z[0], z[1], s) for z in zip(data[cond], data[cont])]
        sizes_processed.loc[len(sizes_processed)] = (prefix + cond, sizes.loc[cond]['size'])
    return scores, sizes_processed


def post_process_callback(path):
    try:
        print('RESULT', path)
        signals_tests.process(path)
        signals_visualize.process(path)
    except Exception as e:
        print(e, file=sys.stderr)


def error_callback(e):
    print(e, file=sys.stderr)


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

    process(data_path, sizes_path, peaks_sizes_path)


if __name__ == "__main__":
    main()
