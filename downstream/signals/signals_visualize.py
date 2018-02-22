import getopt
import math
import sys
from enum import Enum

import matplotlib
import numpy as np
import pandas as pd
from sklearn import preprocessing
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression

from downstream.aging import *

# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # nopep8
import seaborn as sns  # nopep8

# Explicitly setup style
plt.style.use('seaborn-darkgrid')


def pca_separation_fit_error(x_r, y):
    # Fit logistic regression and compute fit error
    lr = LogisticRegression()
    lr.fit(x_r, y)
    p_y = [1 if x[0] < 0.5 else 0 for x in lr.predict_proba(x_r)]
    error = np.sum(np.abs(np.subtract(y, p_y)))

    return error


def pca_signal(signal, scale):
    # Data may contain "Inf" or "NaN" values for some rages, let's just skip
    # such values otherwise PCA will fail
    with pd.option_context('mode.use_inf_as_null', True):
        df = signal.dropna(how="any", axis=0).T
    # Use scaled PCA anyway
    if scale:
        df = preprocessing.scale(df)
    pca = PCA(n_components=2)
    x_r = pca.fit_transform(df)
    return pca, x_r


def signal_pca_plot(signal, title, ax, *,
                    scale=True, groups=None, show_error=True, show_names=True):
    columns = signal.columns
    if groups is None:
        groups = [group(d) for d in columns]
    pca, x_r = pca_signal(signal, scale)

    y = [0 if g == YOUNG else 1 for g in groups]
    if show_error and len(set(y)) > 1:
        error = pca_separation_fit_error(x_r, y)
    else:
        error = None

    for g in set(groups):
        group_filter = np.asarray([g == n for n in groups])
        ax.scatter(x_r[group_filter, 0], x_r[group_filter, 1],
                   color=g.color, alpha=.8, label=g.name)

    for c, g, x, y in zip(columns, groups, x_r[:, 0], x_r[:, 1]):
        ax.annotate(donor(c) if show_names else '',
                    xy=(x, y),
                    xytext=(5, 0),
                    color=g.color,
                    textcoords='offset points',
                    ha='right', va='bottom')
    (var1, var2) = pca.explained_variance_ratio_
    pc1_var = 0.0 if math.isnan(var1) else int(var1 * 100)
    pc2_var = 0.0 if math.isnan(var2) else int(var2 * 100)
    ax.set_xlabel('PC1 {}%'.format(pc1_var))
    ax.set_ylabel('PC2 {}%'.format(pc2_var))
    if error is not None:
        ax.set_title('{} error: {}'.format(title, error))
    else:
        ax.set_title(title)
    return error


class Plot(Enum):
    SCATTER = 1
    HIST = 2
    MA = 3


def mean_regions(df, title, ax, plot_type):
    """Plots for mean values over OD and YD"""
    ods = [c for c in df.columns if is_od(c)]
    yds = [c for c in df.columns if is_yd(c)]

    # Use LOG mean signal, since this is what people are used to
    signal = pd.DataFrame()
    signal["O"] = np.log(df[ods].mean(axis=1)) / np.log(10)
    signal["Y"] = np.log(df[yds].mean(axis=1)) / np.log(10)
    # Fix NA
    signal.loc[~np.isfinite(signal["O"]), "O"] = 0.0
    signal.loc[~np.isfinite(signal["Y"]), "Y"] = 0.0

    if plot_type == Plot.MA:
        signal["M"] = signal["O"] - signal["Y"]
        signal["A"] = 0.5 * (signal["O"] + signal["Y"])
        ax.scatter(signal["A"], signal["M"], alpha=.3, s=1)
        ax.set_xlabel("A")
        ax.set_ylabel("M")

        xmin = np.min(ax.get_xlim())
        xmax = np.max(ax.get_xlim())
        ax.plot([xmin, xmax], [0, 0], c="red", alpha=0.75, lw=1, ls='dotted')
        ax.set_xlim([xmin, xmax])

    elif plot_type == Plot.HIST:
        ax.hist(signal["O"], color=OLD.color, bins=100, alpha=0.3, label="log10(O)")
        ax.hist(signal["Y"], color=YOUNG.color, bins=100, alpha=0.3, label="log10(Y)")
        ax.legend()
    else:
        ax.scatter(signal["O"], signal["Y"], alpha=.3, s=1)
        ax.set_xlabel("log10(O)")
        ax.set_ylabel("log10(Y)")
        # x = y
        lims = [
            np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
            np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
        ]
        # Line x = y
        ax.plot(lims, lims, 'r-', alpha=0.75, lw=1, ls='dotted')
        ax.set_aspect('equal')
        ax.set_xlim(lims)
        ax.set_ylim(lims)

    ax.set_title(title)


def mean_boxplots(df, title, ax):
    """Plot mean values for individual donors"""
    signal = df.mean(axis=1).to_frame('value')
    groups = [group(n) for n in signal.index]
    signal['group'] = [g.name for g in groups]
    all_groups = [g for g in reversed(sorted(set(groups)))]
    all_groups_names = [g.name for g in all_groups]
    sns.boxplot(x='group', y='value', data=signal, palette='Set3',
                linewidth=1.0, order=all_groups_names, ax=ax)
    sns.swarmplot(x='group', y='value', data=signal, color='.25', order=all_groups_names, ax=ax)

    for i, g in enumerate(all_groups):
        group_data = signal[signal['group'] == g.name]
        for j, label in enumerate(group_data.index):
            ax.annotate(donor(label),
                        xy=(i, group_data.iloc[j, :]['value']),
                        xytext=(5, 0),
                        color=g.color,
                        textcoords='offset points')
    ax.set_title(title)


def visualize(f, signal_type):
    print('Visualizing', signal_type, f)
    df = pd.read_csv(f, sep='\t').drop(['chr', 'start', 'end'], axis=1)
    signal = df[[c for c in df.columns.values if not is_input(c)]]
    ods = [c for c in signal.columns if is_od(c)]
    yds = [c for c in signal.columns if is_yd(c)]
    if ods and yds:
        plt.figure(figsize=(36, 6))
        fit_error = signal_pca_plot(signal, title=signal_type, ax=plt.subplot(1, 5, 1))
        mean_regions(df, title='O vs Y {}'.format(signal_type), ax=plt.subplot(1, 5, 2),
                     plot_type=Plot.SCATTER)
        mean_regions(df, title='MA log10 {}'.format(signal_type), ax=plt.subplot(1, 5, 3),
                     plot_type=Plot.MA)
        mean_regions(df, title='Histogram {}'.format(signal_type), ax=plt.subplot(1, 5, 4),
                     plot_type=Plot.HIST)
        mean_boxplots(signal.T, title=signal_type, ax=plt.subplot(1, 5, 5))
        plt.savefig(re.sub('.tsv', '.png', f))
        plt.close()

        # Save pca fit errors to file
        pd.DataFrame(data=[[fit_error]]).to_csv(re.sub('.tsv', '_pca_fit_error.csv', f),
                                                index=None, header=False)
    else:
        plt.figure(figsize=(12, 6))
        signal_pca_plot(signal, title=signal_type, ax=plt.subplot(1, 2, 1))
        mean_boxplots(signal.T, title=signal_type, ax=plt.subplot(1, 2, 2))
        plt.savefig(re.sub('.tsv', '.png', f))
        plt.close()


def process(path):
    visualize(path, re.sub('(.*_)|(\\.tsv$)'.format(id), '', path))


def main():
    argv = sys.argv
    opts, args = getopt.getopt(argv[1:], "h", ["help"])
    if len(args) != 1:
        print("ARGUMENTS:  <data.tsv>\n"
              "<data.tsv> - processed signal file to visualize and build PCA"
              "ARGS: " + ",".join(args))
        sys.exit(1)

    path = args[0]
    process(path)


if __name__ == "__main__":
    main()
