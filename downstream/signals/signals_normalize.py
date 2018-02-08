#!/usr/bin/env python
# Author: oleg.shpynov@jetbrains.com

import math
import multiprocessing
import tempfile

import numpy as np
import pandas as pd
from sklearn import preprocessing

from downstream.aging import *
from downstream.signals import signals_visualize
from scripts.util import *


def process(data_path, sizes_path, peaks_sizes_path, *, processes=4):
    pool = multiprocessing.Pool(processes=processes)
    pool.apply_async(raw_normalization,
                     args=(data_path,),
                     error_callback=error_callback)
    pool.apply_async(rpm_normalization,
                     args=(data_path, sizes_path),
                     error_callback=error_callback)

    pool.apply_async(diffbind_tmm_minus_full,
                     args=(data_path, sizes_path),
                     error_callback=error_callback)
    pool.apply_async(diffbind_tmm_reads_full_cpm,
                     args=(data_path, sizes_path),
                     error_callback=error_callback)

    if peaks_sizes_path and os.path.exists(peaks_sizes_path):
        pool.apply_async(frip_normalization,
                         args=(data_path, sizes_path, peaks_sizes_path),
                         error_callback=error_callback)
        pool.apply_async(diffbind_tmm_reads_effective_cpm,
                         args=(data_path, peaks_sizes_path),
                         error_callback=error_callback)
    pool.close()
    pool.join()


def processing_chipseq(loaded):
    return not [n for n in loaded['name'] if re.match('.*(meth|trans|mirna).*', n)]


def raw_normalization(data_path):
    """Raw signal with standard scaling"""
    raw_path = re.sub('.tsv', '_raw.tsv', data_path)
    z_path = re.sub('.tsv', '_rawz.tsv', data_path)
    if os.path.exists(raw_path) and os.path.exists(z_path):
        return

    loaded = pd.read_csv(data_path, sep='\t',
                         names=('chr', 'start', 'end', 'coverage', 'mean0', 'mean', 'name'))
    print('Processing RAW signal')
    raw_data = pd.pivot_table(loaded, index=['chr', 'start', 'end'],
                              columns='name',
                              values='coverage' if processing_chipseq(loaded) else 'mean',
                              fill_value=0)
    if not os.path.exists(raw_path):
        raw_data.to_csv(raw_path, sep='\t')
        print('Saved RAW signal to {}'.format(raw_path))
        signals_visualize.process(raw_path)

    if not os.path.exists(z_path):
        print("Processing Z normalization")
        process_z(z_path, raw_data)
        signals_visualize.process(z_path)


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
        signals_visualize.process(rpm_path)

    if not os.path.exists(rpkm_path):
        print('Processing RPKM normalization')
        data['rpk'] = (data['end'] - data['start']) / 1000.0
        data['rpkm'] = data['rpm'] / data['rpk']
        save_signal(rpkm_path, data, 'rpkm', 'Saved RPKM')
        signals_visualize.process(rpkm_path)


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


def frip_normalization(data_path, sizes_path, peaks_sizes_path):
    """Normalization on library depths, FRIPs"""
    frip_pm_path = re.sub('.tsv', '_fripm.tsv', data_path)
    fripz_path = re.sub('.tsv', '_fripz.tsv', data_path)
    if os.path.exists(frip_pm_path) and os.path.exists(fripz_path):
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

    # Subtract number number of noisy reads
    frip_df = pd.DataFrame()
    for column in not_input_columns:
        frip_df[column] = np.maximum(0, data[column] - input_coef[column] * data[input_column])

    # Per million reads normalization
    if not os.path.exists(frip_pm_path):
        frip_pm = pd.DataFrame()
        for column in not_input_columns:
            frip_pm[column] = frip_df[column] * 1000000.0 / counts['data'][column]
        frip_pm[not_input_columns].to_csv(frip_pm_path, sep='\t')
        print('{} to {}'.format('Saved FRIP per million reads normalized', frip_pm_path))
        signals_visualize.process(frip_pm_path)

    # Z normalization
    if not os.path.exists(fripz_path):
        process_z(fripz_path, frip_df)
        signals_visualize.process(fripz_path)


def diffbind_tmm_reads_effective_cpm(data_path, peaks_sizes_path):
    scores_tmm_reads_effective_cpm_path = \
        re.sub('.tsv', '_diffbind_tmm_reads_effective_cpm.tsv', data_path)
    if os.path.exists(scores_tmm_reads_effective_cpm_path):
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
    scores_tmm_reads_effective_cpm = process_tmm(data, peaks_sizes) * 1000000.0
    scores_tmm_reads_effective_cpm.to_csv(scores_tmm_reads_effective_cpm_path, sep='\t')
    print('Saved Diffbind DBA_SCORE_TMM_READS_EFFECTIVE_CPM to {}'.format(scores_tmm_reads_effective_cpm_path))
    signals_visualize.process(scores_tmm_reads_effective_cpm_path)


def diffbind_tmm_reads_full_cpm(data_path, sizes_path):
    scores_tmm_reads_full_cpm_path = \
        re.sub('.tsv', '_diffbind_tmm_reads_full_cpm.tsv', data_path)
    if os.path.exists(scores_tmm_reads_full_cpm_path):
        return

    loaded = pd.read_csv(data_path, sep='\t',
                         names=('chr', 'start', 'end', 'coverage', 'mean0', 'mean', 'name'))
    if not processing_chipseq(loaded):
        return

    print('Processing DBA_SCORE_TMM_READS_FULL_CPM')
    sizes = pd.read_csv(sizes_path, sep='\t', names=('name', 'size'), index_col='name')
    data = pd.pivot_table(loaded, index=['chr', 'start', 'end'],
                          columns='name', values='coverage', fill_value=0)
    scores_tmm_reads_full_cpm = process_tmm(data, sizes) * 1000000.0
    scores_tmm_reads_full_cpm.to_csv(scores_tmm_reads_full_cpm_path, sep='\t')
    print('Saved Diffbind DBA_SCORE_TMM_READS_FULL_CPM to {}'.format(scores_tmm_reads_full_cpm_path))
    signals_visualize.process(scores_tmm_reads_full_cpm_path)


def diffbind_tmm_minus_full(data_path, sizes_path):
    scores_tmm_minus_full_path = \
        re.sub('.tsv', '_diffbind_tmm_minus_full.tsv', data_path)
    if os.path.exists(scores_tmm_minus_full_path):
        return

    loaded = pd.read_csv(data_path, sep='\t',
                         names=('chr', 'start', 'end', 'coverage', 'mean0', 'mean', 'name'))
    if not processing_chipseq(loaded):
        return
    print('Processing DBA_SCORE_TMM_MINUS_FULL')

    sizes = pd.read_csv(sizes_path, sep='\t', names=('name', 'size'), index_col='name')
    data = pd.pivot_table(loaded, index=['chr', 'start', 'end'],
                          columns='name', values='coverage', fill_value=0)
    pairs = [(c, find_input(c, data.columns)) for c in data.columns if not is_input(c)]
    scores_minus_scaled_control = diffbind_scores_minus(data, sizes, pairs)
    scores_tmm_minus_full = \
        process_tmm(scores_minus_scaled_control, sizes) * float(np.mean(sizes))
    scores_tmm_minus_full.to_csv(scores_tmm_minus_full_path, sep='\t')
    print('Saved Diffbind DBA_SCORE_TMM_MINUS_FULL to {}'.format(scores_tmm_minus_full_path))
    signals_visualize.process(scores_tmm_minus_full_path)


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
        scores[condition] = [score(z[0], z[1], scale) for z in zip(data[condition], data[control])]
    scores.index = data.index
    return scores


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
