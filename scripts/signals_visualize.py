import matplotlib

# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

import getopt
import math
import os
import re
import subprocess
import sys
import tempfile
from collections import namedtuple
from enum import Enum
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn import preprocessing
from sklearn.decomposition import PCA


class Normalization(Enum):
    NONE = 1
    SCALED = 2


Group = namedtuple('Group', 'name color prefix')
OD_GROUP = Group('OD', 'blue', '')
YD_GROUP = Group('YD', 'red', '')


def signal_pca(x0,
               title,
               groups=None,
               scale=Normalization.NONE):
    if groups is None:
        groups = [OD_GROUP if re.match('.*od\\d+.*', r, flags=re.IGNORECASE) else YD_GROUP for r in x0.index]
    x = x0
    if scale == Normalization.SCALED:
        x = preprocessing.scale(x0)

    pca = PCA(n_components=2)
    x_r = pca.fit_transform(x)
    for g in set(groups):
        group_filter = np.asarray([g == n for n in groups])
        plt.scatter(x_r[group_filter, 0], x_r[group_filter, 1], color=g.color, alpha=.8, label=g.name)

    for g, label, x, y in zip(groups,
                              [re.search('[yo]d\\d+', n, flags=re.IGNORECASE).group(0) for n in x0.index],
                              x_r[:, 0], x_r[:, 1]):
        plt.annotate(g.prefix + label,
                     xy=(x, y),
                     xytext=(5, 0),
                     color=g.color,
                     textcoords='offset points',
                     ha='right', va='bottom')
    plt.title(title)
    plt.xlabel('PC1 {}%'.format(int(pca.explained_variance_ratio_[0] * 100)))
    plt.ylabel('PC2 {}%'.format(int(pca.explained_variance_ratio_[1] * 100)))


def signal_pca_all(x, title, groups=None):
    """Plot all the scaled variants of PCA"""
    plt.figure(figsize=(20, 5))
    plt.subplot(1, 4, 1)
    signal_pca(x, title, groups=groups)
    plt.subplot(1, 4, 2)
    signal_pca(x, 'SCALED {}'.format(title), scale=Normalization.SCALED, groups=groups)
    plt.subplot(1, 4, 3)
    signal_pca(np.log1p(x), 'LOG {}'.format(title), groups=groups)
    plt.subplot(1, 4, 4)
    signal_pca(np.log1p(x), 'SCALED LOG {}'.format(title), scale=Normalization.SCALED, groups=groups)


class Plot(Enum):
    SCATTER = 1
    HISTOGRAM = 2
    MA = 3


def mean_regions(df, title, ax, plot_type):
    """Plots for mean values over OD and YD"""
    ods = [c for c in df.columns if
           re.match('.*od\\d+.*', c, flags=re.IGNORECASE) and not re.match('.*input.*', c, flags=re.IGNORECASE)]
    yds = [c for c in df.columns if
           re.match('.*yd\\d+.*', c, flags=re.IGNORECASE) and not re.match('.*input.*', c, flags=re.IGNORECASE)]

    signal = pd.DataFrame()
    signal["ODS"] = df[ods].mean(axis=1)
    signal["YDS"] = df[yds].mean(axis=1)

    if plot_type == Plot.MA:
        signal["M"] = np.log1p(signal["ODS"]) - np.log1p(signal["YDS"])
        signal["A"] = 0.5 * (np.log1p(signal["ODS"]) + np.log1p(signal["YDS"]))
        ax.scatter(signal["A"], signal["M"], alpha=.3, s=1)
        ax.set_xlabel("A")
        ax.set_ylabel("M")

        xmin = np.min(ax.get_xlim())
        xmax = np.max(ax.get_xlim())
        ax.plot([xmin, xmax], [0, 0], c="red", alpha=0.75, lw=1, ls='dotted')
        ax.set_xlim([xmin, xmax])

    elif plot_type == Plot.HISTOGRAM:
        signal["ODS"] = np.log1p(signal["ODS"]) / np.log(10)
        signal["YDS"] = np.log1p(signal["YDS"]) / np.log(10)

        ax.hist(signal["ODS"], color=OD_GROUP.color, bins=100, alpha=0.3, label="ODS")
        ax.hist(signal["YDS"], color=YD_GROUP.color, bins=100, alpha=0.3, label="YDS")
        ax.legend()
    else:
        ax.scatter(signal["ODS"], signal["YDS"], alpha=.3, s=1)
        ax.set_xlabel("mean ODS")
        ax.set_ylabel("mean YDS")
        # x = y
        lims = [
            np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
            np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
        ]
        # now plot both limits against eachother
        ax.plot(lims, lims, 'r-', alpha=0.75, lw=1, ls='dotted')
        ax.set_aspect('equal')
        ax.set_xlim(lims)
        ax.set_ylim(lims)

    ax.set_title(title)


def mean_boxplots(df, title, ax):
    """Plot mean values for individual donors"""
    signal = df.mean(axis=1).to_frame("value")
    signal.index = [re.search('[yo]d\\d+', n, flags=re.IGNORECASE).group(0) for n in signal.index]
    signal["age"] = "ODS"
    signal.loc[signal.index.str.startswith("y"), "age"] = "YDS"

    age_labels = list(reversed(sorted(list(set(signal['age'])))))
    sns.boxplot(x="age", y="value", data=signal, palette="Set3", linewidth=1.0, order=age_labels, ax=ax)
    sns.swarmplot(x="age", y="value", data=signal, color=".25", order=age_labels, ax=ax)

    for i, age_label in enumerate(age_labels):
        age_data = signal[signal['age'] == age_label]
        for j, label in enumerate(age_data.index):
            ax.annotate(label, xy=(i, age_data.iloc[j, :]['value']),
                        xytext=(5, 0),
                        color=OD_GROUP.color if age_label == "YDS" else YD_GROUP.color,
                        textcoords='offset points')

    ax.set_title(title)
    return signal


# See documents on how to compute scores
# https://docs.google.com/document/d/1zH5cw5Zal546xkoFFCVqhhYmf3742efhddz5cqpD9PQ/edit?usp=sharing
def score(cond, cont, scale):
    if scale > 1:
        scale = 1
    if scale != 0:
        cont = math.ceil(cont * scale)
    # According to pv.get_reads() function, if reads number is < 1 than it should be 1
    return max(1, cond - cont)


def process_scores(df, sizes, records):
    scores_processed = pd.DataFrame()
    sizes_processed = pd.DataFrame(columns=['name', 'size'])
    for cond, cont, g in records:
        scale = sizes.loc[cond]['size'] / sizes.loc[cont]['size']
        if g is None:
            prefix = ''
        else:
            prefix = g.prefix
        scores_processed[prefix + cond] = [score(z[0], z[1], scale) for z in zip(df[cond], df[cont])]
        sizes_processed.loc[len(sizes_processed)] = (prefix + cond, sizes.loc[cond]['size'])
    return scores_processed, sizes_processed


TMM_R_PATH = os.path.dirname(os.path.realpath(__file__)) + '/../R/tmm.R'


def process(folder, id):
    print('Processing signal data')
    file = folder + '/{0}/{0}_{1}.tsv'
    for normalization in ['raw', 'rpm', 'rpkm', 'rpm_peaks']:
        f = file.format(id, normalization)
        try:
            print(f)
            df = pd.read_csv(f, sep='\t')
            od_input = [c for c in df.columns.values if re.match('.*od.*input.*', c, flags=re.IGNORECASE)][0]
            yd_input = [c for c in df.columns.values if re.match('.*yd.*input.*', c, flags=re.IGNORECASE)][0]
            signal = df.drop(['chr', 'start', 'end', od_input, yd_input], axis=1)
            plt.figure(figsize=(20, 5))
            mean_regions(df, title=normalization, ax=plt.subplot(1, 4, 1), plot_type=Plot.SCATTER)
            mean_regions(df, title='MA {}'.format(normalization), ax=plt.subplot(1, 4, 2), plot_type=Plot.MA)
            mean_regions(df, title='LOG {}'.format(normalization), ax=plt.subplot(1, 4, 3), plot_type=Plot.HISTOGRAM)
            means = mean_boxplots(signal.T, title=normalization, ax=plt.subplot(1, 4, 4))
            plt.savefig(re.sub('.tsv', '.png', f))
            plt.close()

            # Save means signal to df
            pd.DataFrame(means['value']).to_csv(re.sub('.tsv', '_data.csv', f), index=True, header=None)

            signal_pca_all(signal.T, normalization)
            plt.savefig(re.sub('.tsv', '_pca.png', f))
            plt.close()

        except FileNotFoundError:
            print('File not found: {}'.format(f))

    print('Processing diffbind normalization')
    f = folder + '/{0}/{0}_raw.tsv'.format(id)
    try:
        df = pd.read_csv(f, sep='\t')
        od_input = [c for c in df.columns.values if re.match('.*od.*input.*', c, flags=re.IGNORECASE)][0]
        yd_input = [c for c in df.columns.values if re.match('.*yd.*input.*', c, flags=re.IGNORECASE)][0]
        ods = [c for c in df.columns if
               re.match('.*od\\d+.*', c, flags=re.IGNORECASE) and not re.match('.*input.*', c, flags=re.IGNORECASE)]
        yds = [c for c in df.columns if
               re.match('.*yd\\d+.*', c, flags=re.IGNORECASE) and not re.match('.*input.*', c, flags=re.IGNORECASE)]

        sizes = pd.read_csv(folder + '/{0}/sizes.tsv'.format(id), sep='\t', names=('name', 'size'))
        sizes.index = sizes['name']
        sizes = sizes.drop('name', axis=1)
        sizes['size'] = sizes['size'] / 1000000

        records = [(d, od_input, OD_GROUP) for d in ods] + [(d, yd_input, YD_GROUP) for d in yds]
        scores, lib_sizes = process_scores(df, sizes, records)

        plt.figure(figsize=(20, 5))
        mean_regions(scores, title='diffbind score', ax=plt.subplot(1, 4, 1), plot_type=Plot.SCATTER)
        mean_regions(scores, title='MA diffbind score', ax=plt.subplot(1, 4, 2), plot_type=Plot.MA)
        mean_regions(scores, title='LOG diffbind score', ax=plt.subplot(1, 4, 3), plot_type=Plot.HISTOGRAM)
        means = mean_boxplots(scores.T, title='diffbind', ax=plt.subplot(1, 4, 4))
        plt.savefig(re.sub('_raw.tsv', '_scores.png', f))
        plt.close()

        # Save means signal to df
        pd.DataFrame(means['value']).to_csv(re.sub('.tsv', '_scores_data.csv', f), index=True, header=None)

        groups = [p[2] for p in records]
        signal_pca_all(scores.T, 'Scores', groups=groups)
        plt.savefig(re.sub('_raw.tsv', '_scores_pca.png', f))
        plt.close()

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

        # counts <- counts * mean(res$samples$lib.size)
        scores_tmm_full = scores_tmm.T * sizes.loc[ods + yds]['size'].mean()
        plt.figure(figsize=(20, 5))
        plt.subplot(1, 4, 1)
        signal_pca(scores_tmm.T, 'TMM DBA_SCORE_TMM_MINUS_FULL_CPM', groups=groups)
        plt.subplot(1, 4, 2)
        signal_pca(np.log1p(scores_tmm.T), 'LOG TMM DBA_SCORE_TMM_MINUS_FULL_CPM', groups=groups)
        plt.subplot(1, 4, 3)
        signal_pca(scores_tmm_full, 'TMM DBA_SCORE_TMM_MINUS_FULL', groups=groups)
        plt.subplot(1, 4, 4)
        signal_pca(np.log1p(scores_tmm_full), 'LOG TMM DBA_SCORE_TMM_MINUS_FULL', groups=groups)
        plt.savefig(re.sub('_raw.tsv', '_scores_tmm_pca.png', f))
        plt.close()

    except FileNotFoundError:
        print('File not found: {}'.format(f))


help_message = 'ARGUMENTS: <folder with BW or BAM> <id>'


def usage():
    print(help_message)


def main():
    argv = sys.argv
    opts, args = getopt.getopt(argv[1:], "h", ["help"])
    if len(args) != 2:
        usage()
        sys.exit(1)

    folder = args[0]
    id = args[1]

    print('Processing signal_visualize.py {} {}'.format(folder, id))
    process(folder, id)


if __name__ == "__main__":
    main()
