#!/usr/bin/env python
# Author: oleg.shpynov@jetbrains.com
import sys
import math
import multiprocessing
import tempfile

import numpy as np
import pandas as pd
from sklearn import linear_model

from downstream.signals import signals_visualize
from scripts.util import *


def process(data_path, sizes_path, peaks_sizes_path, *, processes=7):
    pool = multiprocessing.Pool(processes=processes)
    pool.apply_async(raw_normalization,
                     args=(data_path,),
                     error_callback=error_callback)
    pool.apply_async(rpm_normalization,
                     args=(data_path, sizes_path),
                     error_callback=error_callback)

    pool.apply_async(ma_normalization,
                     args=(data_path,),
                     error_callback=error_callback)

    pool.apply_async(diffbind_tmm_minus_full,
                     args=(data_path, sizes_path),
                     error_callback=error_callback)
    pool.apply_async(diffbind_tmm_reads_full_cpm,
                     args=(data_path, sizes_path),
                     error_callback=error_callback)

    if peaks_sizes_path:
        if not os.path.exists(peaks_sizes_path):
            print("File not found:", peaks_sizes_path, file=sys.stderr)
        else:
            pool.apply_async(frip_normalization,
                             args=(data_path, sizes_path, peaks_sizes_path),
                             error_callback=error_callback)
            pool.apply_async(diffbind_tmm_reads_effective_cpm,
                             args=(data_path, peaks_sizes_path),
                             error_callback=error_callback)
    pool.close()
    pool.join()


def processing_chipseq(loaded):
    return not any(re.match('.*(meth|trans|mirna).*', n) for n in loaded['name'])


def raw_normalization(data_path):
    """Raw signal with Q scaling"""
    raw_path = re.sub('.tsv', '_raw.tsv', data_path)
    q_path = re.sub('.tsv', '_rawq.tsv', data_path)
    q_no_ref_path = re.sub('.tsv', '_rawqnr.tsv', data_path)
    if os.path.exists(raw_path) and os.path.exists(q_path) and os.path.exists(q_no_ref_path):
        return

    loaded = pd.read_csv(data_path, sep='\t',
                         names=('chr', 'start', 'end', 'coverage', 'mean0', 'mean', 'name'))
    print('Processing RAW signal')
    raw_data = pd.pivot_table(loaded, index=['chr', 'start', 'end'],
                              columns='name',
                              values='coverage' if processing_chipseq(loaded) else 'mean',
                              fill_value=0)

    raw_data.to_csv(raw_path, sep='\t')
    print('Saved RAW signal to {}'.format(raw_path))
    visualize_safe(raw_path)

    print("Processing RAW Q normalization")
    process_quantile(q_path, raw_data)
    visualize_safe(q_path)

    print("Processing RAW Q_NO_REF normalization")
    process_quantile_no_ref(q_no_ref_path, raw_data)
    visualize_safe(q_no_ref_path)


def rpm_normalization(data_path, sizes_path):
    """RPM/RPKM normalization"""
    rpm_path = re.sub('.tsv', '_rpm.tsv', data_path)
    rpkm_path = re.sub('.tsv', '_rpkm.tsv', data_path)
    if os.path.exists(rpm_path) and os.path.exists(rpkm_path):
        return

    loaded = pd.read_csv(data_path, sep='\t',
                         names=('chr', 'start', 'end', 'coverage', 'mean0', 'mean', 'name'))
    if not processing_chipseq(loaded):
        return
    sizes = pd.read_csv(sizes_path, sep='\t', names=('name', 'size'), index_col='name')
    sizes_pm = sizes / 1000000.0
    sizes_pm.columns = ['size_pm']
    data = pd.merge(loaded, sizes_pm, left_on="name", how='left', right_index=True)

    if not os.path.exists(rpm_path):
        print('Processing RPM normalization')
        data['rpm'] = data['coverage'] / data['size_pm']
        save_signal(rpm_path, data, 'rpm', 'Saved RPM')
        visualize_safe(rpm_path)

    if not os.path.exists(rpkm_path):
        print('Processing RPKM normalization')
        data['rpk'] = (data['end'] - data['start']) / 1000.0
        data['rpkm'] = data['rpm'] / data['rpk']
        save_signal(rpkm_path, data, 'rpkm', 'Saved RPKM')
        visualize_safe(rpkm_path)


def save_signal(path, data, signal_type, msg):
    pivot_df = pd.pivot_table(data, index=['chr', 'start', 'end'],
                              columns='name', values=signal_type, fill_value=0)
    pivot_df.to_csv(path, sep='\t')
    print('{} to {}'.format(msg, path))


def process_quantile(output, data):
    # Normalize everything to the first track
    signal_columns = [c for c in data.columns if not is_input(c)]

    target_column = signal_columns[0]
    print("Quantile normalization to", target_column)
    target = data[target_column]
    target_sorted = np.sort(target)
    signal_df = data[signal_columns]
    ranks_df = signal_df.rank(method='first').astype(int)
    quantile_df = pd.DataFrame.from_dict(
        {c: target_sorted[ranks_df[c] - 1] for c in signal_columns}
    )
    quantile_df.index = data.index
    quantile_df.to_csv(output, sep='\t')
    print('{} to {}'.format('Saved QUANTILE', output))


def process_quantile_no_ref(output, data):
    # To quantile normalize two or more distributions to each other, without a reference
    # distribution, sort as before, then set to the average (usually, arithmetic mean) of
    #  the distributions. So the highest value in all cases becomes the mean of the highest
    # values, the second highest value becomes the mean of the second highest values, and so on.
    #
    # See https://en.wikipedia.org/wiki/Quantile_normalization
    signal_columns = [c for c in data.columns if not is_input(c)]

    signal_df = data[signal_columns]
    ranks_mean = signal_df.stack().groupby(signal_df.rank(method='first').astype(int).stack())\
        .mean()
    quantile_df = signal_df.rank(method='min').astype(int).stack().map(ranks_mean).unstack()
    quantile_df.to_csv(output, sep='\t')
    print('Saved QUANTILE_NO_REF to {}'.format(output))


def ma_normalization(data_path):
    """Normalization based on MA plot"""
    manm_path = re.sub('.tsv', '_manorm.tsv', data_path)

    if os.path.exists(manm_path):
        return

    loaded = pd.read_csv(data_path, sep='\t',
                         names=('chr', 'start', 'end', 'coverage', 'mean0', 'mean', 'name'))

    if not processing_chipseq(loaded):
        return

    print('Processing MA normalization')
    data = pd.pivot_table(loaded, index=['chr', 'start', 'end'],
                          columns='name', values='coverage', fill_value=0)
    not_input_columns = [c for c in data.columns if not is_input(c)]

    mean = data[not_input_columns].values.mean(axis=1)

    ma_df = pd.DataFrame()
    for column in not_input_columns:
        values = data[column].copy()

        indexes = np.logical_and(values > 0, mean > 0).values
        x = values[indexes].values
        y = mean[indexes]

        log_x = np.log(x)
        log_y = np.log(y)

        m = log_x - log_y
        a = 0.5 * (log_x + log_y)

        regr = linear_model.LinearRegression()
        a_rs = a.reshape((len(a), 1))
        regr.fit(a_rs, m)
        m_pred = regr.predict(a_rs)

        x_new = np.exp(log_x - m_pred)

        values[indexes] = x_new

        ma_df[column] = values

    ma_df.to_csv(manm_path, sep='\t')
    print('{} to {}'.format('Saved MA normalized', manm_path))
    visualize_safe(manm_path)


def frip_normalization(data_path, sizes_path, peaks_sizes_path):
    """Normalization on library depths, FRIPs"""
    frip_pm_path = re.sub('.tsv', '_fripm.tsv', data_path)
    frip_pm_1kbp_path = re.sub('.tsv', '_fripm_1kbp.tsv', data_path)
    if os.path.exists(frip_pm_path) and os.path.exists(frip_pm_1kbp_path):
        return

    loaded = pd.read_csv(data_path, sep='\t',
                         names=('chr', 'start', 'end', 'coverage', 'mean0', 'mean', 'name'))
    if not processing_chipseq(loaded):
        return
    print('Processing FRIP normalization')
    sizes = pd.read_csv(sizes_path, sep='\t', names=('name', 'size'), index_col='name')
    data = pd.pivot_table(loaded, index=['chr', 'start', 'end'],
                          columns='name', values='coverage', fill_value=0)
    peaks_sizes = pd.read_csv(peaks_sizes_path, sep='\t', names=('name', 'tags_in_peaks'),
                              index_col='name')
    counts = pd.merge(sizes, peaks_sizes, left_index=True, right_index=True)

    # Here we assume that we have single input for everything
    input_column = [c for c in data.columns if is_input(c)][0]
    not_input_columns = [c for c in data.columns if not is_input(c)]
    counts['input_tags'] = counts.loc[input_column, 'size']
    counts['input_tags_in_peaks'] = counts.loc[input_column, 'tags_in_peaks']

    # Omit input columns
    counts = counts[~ counts.index.str.contains("input")]

    no_peaks_coverages = counts['size'] - counts['tags_in_peaks']
    no_peaks_input_coverages = counts['input_tags'] - counts['input_tags_in_peaks']

    # Input coefficient shows the ratio of signal track noise to input track level.
    # The smaller coefficient the better ChIP-Seq tracks we deal with.
    input_coef = no_peaks_coverages / no_peaks_input_coverages
    counts['input_coef'] = input_coef
    counts['data'] = counts['tags_in_peaks'] - input_coef * counts['input_tags_in_peaks']

    # Subtract number of noisy reads
    frip_df = pd.DataFrame()
    for column in not_input_columns:
        frip_df[column] = np.maximum(0, data[column] - input_coef[column] * data[input_column])

    # Per million reads normalization
    frip_pm = pd.DataFrame()
    for column in not_input_columns:
        frip_pm[column] = frip_df[column] * 1000000.0 / counts['data'][column]
    frip_pm[not_input_columns].to_csv(frip_pm_path, sep='\t')
    print('{} to {}'.format('Saved FRIP per million reads normalized', frip_pm_path))
    visualize_safe(frip_pm_path)

    # Per 1kbp normalization
    frip_pm_1kbp = pd.DataFrame(frip_pm.to_records())
    for column in not_input_columns:
        frip_pm_1kbp[column] = \
            frip_pm_1kbp[column] * 1000.0 / (frip_pm_1kbp['end'] - frip_pm_1kbp['start'])
    frip_pm_1kbp.index = frip_pm.index

    frip_pm_1kbp[not_input_columns].to_csv(frip_pm_1kbp_path, sep='\t')
    print('{} to {}'.format('Saved FRIP per 1kbp per million reads normalized', frip_pm_1kbp_path))
    visualize_safe(frip_pm_1kbp_path)


def diffbind_tmm_reads_effective_cpm(data_path, peaks_sizes_path):
    path = re.sub('.tsv', '_diffbind_tmm_reads_effective_cpm.tsv', data_path)
    path_1kbp = re.sub('.tsv', '_diffbind_tmm_reads_effective_cpm_1kbp.tsv', data_path)
    if os.path.exists(path) and os.path.exists(path_1kbp):
        return

    loaded = pd.read_csv(data_path, sep='\t',
                         names=('chr', 'start', 'end', 'coverage', 'mean0', 'mean', 'name'))
    if not processing_chipseq(loaded):
        return

    print('Processing DBA_SCORE_TMM_READS_EFFECTIVE_CPM')
    data = pd.pivot_table(loaded, index=['chr', 'start', 'end'],
                          columns='name', values='coverage', fill_value=0)
    peaks_sizes = pd.read_csv(peaks_sizes_path, sep='\t', names=('name', 'tags_in_peaks'),
                              index_col='name')
    scores = process_tmm(data, peaks_sizes) * 1000000.0
    scores.to_csv(path, sep='\t')
    print('Saved Diffbind DBA_SCORE_TMM_READS_EFFECTIVE_CPM to {}'.format(path))
    visualize_safe(path)

    # Per 1kbp normalization
    normalize_1kbp(path_1kbp, scores)
    print('Saved Diffbind DBA_SCORE_TMM_READS_EFFECTIVE_CPM 1kbp to {}'.format(path_1kbp))
    visualize_safe(path_1kbp)


def diffbind_tmm_reads_full_cpm(data_path, sizes_path):
    path = re.sub('.tsv', '_diffbind_tmm_reads_full_cpm.tsv', data_path)
    path_1kbp = re.sub('.tsv', '_diffbind_tmm_reads_full_cpm_1kbp.tsv', data_path)
    if os.path.exists(path) and os.path.exists(path_1kbp):
        return

    loaded = pd.read_csv(data_path, sep='\t',
                         names=('chr', 'start', 'end', 'coverage', 'mean0', 'mean', 'name'))
    if not processing_chipseq(loaded):
        return

    print('Processing DBA_SCORE_TMM_READS_FULL_CPM')
    sizes = pd.read_csv(sizes_path, sep='\t', names=('name', 'size'), index_col='name')
    data = pd.pivot_table(loaded, index=['chr', 'start', 'end'],
                          columns='name', values='coverage', fill_value=0)
    scores = process_tmm(data, sizes) * 1000000.0
    scores.to_csv(path, sep='\t')
    print('Saved Diffbind DBA_SCORE_TMM_READS_FULL_CPM to {}'.format(path))
    visualize_safe(path)

    # Per 1kbp normalization
    normalize_1kbp(path_1kbp, scores)
    print('Saved Diffbind DBA_SCORE_TMM_READS_FULL_CPM 1kbp to {}'.format(path_1kbp))
    visualize_safe(path_1kbp)


def diffbind_tmm_minus_full(data_path, sizes_path):
    path = re.sub('.tsv', '_diffbind_tmm_minus_full.tsv', data_path)
    path_1kbp = re.sub('.tsv', '_diffbind_tmm_minus_full_1kbp.tsv', data_path)
    if os.path.exists(path) and os.path.exists(path_1kbp):
        return

    loaded = pd.read_csv(data_path, sep='\t',
                         names=('chr', 'start', 'end', 'coverage', 'mean0', 'mean', 'name'))
    if not processing_chipseq(loaded):
        return
    print('Processing DBA_SCORE_TMM_MINUS_FULL')

    sizes = pd.read_csv(sizes_path, sep='\t', names=('name', 'size'), index_col='name')
    data = pd.pivot_table(loaded, index=['chr', 'start', 'end'],
                          columns='name', values='coverage', fill_value=0)
    pairs = [(c, find_input_name(c, data.columns))
             for c in data.columns if not is_input(c)]
    scores_minus_scaled_control = diffbind_scores_minus(data, sizes, pairs)
    scores = process_tmm(scores_minus_scaled_control, sizes) * float(np.mean(sizes))
    scores.to_csv(path, sep='\t')
    print('Saved Diffbind DBA_SCORE_TMM_MINUS_FULL to {}'.format(path))
    visualize_safe(path)

    # Per 1kbp normalization
    normalize_1kbp(path_1kbp, scores)
    print('Saved Diffbind DBA_SCORE_TMM_MINUS_FULL 1kbp to {}'.format(path_1kbp))
    visualize_safe(path_1kbp)


def normalize_1kbp(path_1kbp, scores):
    not_input_columns = [c for c in scores.columns if not is_input(c)]
    scores_1kbp = pd.DataFrame(scores[not_input_columns].to_records())
    for column in not_input_columns:
        scores_1kbp[column] = scores_1kbp[column] * 1000.0 /\
                              (scores_1kbp['end'] - scores_1kbp['start'])
    scores_1kbp.to_csv(path_1kbp, sep='\t', index=False)


def process_tmm(data, sizes):
    print('Scores TMM normalization')
    scores_tmpfile = tempfile.NamedTemporaryFile(prefix='scores', suffix='.tsv').name
    data.to_csv(scores_tmpfile, index=False, sep='\t')
    tmm_file = scores_tmpfile.replace('.tsv', '_tmm.tsv')
    sizes_tmpfile = tempfile.NamedTemporaryFile(prefix='sizes', suffix='.tsv').name
    sizes.to_csv(sizes_tmpfile, sep='\t', header=None)
    tmm_R = os.path.dirname(os.path.realpath(__file__)) + '/tmm.R'
    cmd = "Rscript " + tmm_R + " " + scores_tmpfile + " " + sizes_tmpfile + " " + tmm_file
    print(cmd)
    subprocess.run(cmd, shell=True)
    result = pd.read_csv(tmm_file, sep='\t')
    result.index = data.index
    return result


def diffbind_scores_minus(data, sizes, pairs):
    """
    Computes DiffBind score DBA_SCORE_TMM_MINUS_FULL
    See documents on how to compute scores
    https://docs.google.com/document/d/1zH5cw5Zal546xkoFFCVqhhYmf3742efhddz5cqpD9PQ/edit?usp=sharing
    """

    def score(cond, cont, s):
        if s > 1:
            s = 1
        if s != 0:
            cont = math.ceil(cont * s)
        # According to DiffBind pv.get_reads() function (utils.R:239)
        # if reads number is < 1 than it should be 1
        return max(1, cond - cont)

    scores = pd.DataFrame()
    for condition, control in pairs:
        scale = sizes.loc[condition]['size'] / sizes.loc[control]['size']
        scores[condition] = [score(z[0], z[1], scale)
                             for z in zip(data[condition], data[control])]
    scores.index = data.index
    return scores


def visualize_safe(path):
    try:
        signals_visualize.process(path)
    except:
        print("Failed to visualize", path)


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
