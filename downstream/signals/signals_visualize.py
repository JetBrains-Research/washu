import getopt
import matplotlib
from enum import Enum
import math
import numpy as np
import pandas as pd
import sys

from downstream.aging import *
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression

# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # nopep8
import seaborn as sns  # nopep8

# Explicitly setup style
plt.style.use('seaborn-darkgrid')


def signal_pca(x, title, ax):
    groups = [OLD if is_od(r) else YOUNG for r in x.index]
    pca = PCA(n_components=2)
    x_r = pca.fit_transform(x)
    for g in set(groups):
        group_filter = np.asarray([g == n for n in groups])
        ax.scatter(x_r[group_filter, 0], x_r[group_filter, 1],
                   color=g.color, alpha=.8, label=g.name)

    for g, label, x, y in zip(groups, [age(n) for n in x.index], x_r[:, 0], x_r[:, 1]):
        ax.annotate(g.prefix + label,
                    xy=(x, y),
                    xytext=(5, 0),
                    color=g.color,
                    textcoords='offset points',
                    ha='right', va='bottom')
    (var1, var2) = pca.explained_variance_ratio_
    var2 = pca.explained_variance_ratio_[0]
    pc1_var = 0.0 if math.isnan(var1) else int(var1 * 100)
    pc2_var = 0.0 if math.isnan(var2) else int(var2 * 100)
    ax.set_xlabel('PC1 {}%'.format(pc1_var))
    ax.set_ylabel('PC2 {}%'.format(pc2_var))

    # Fit logistic regression and compute fit error
    y = [0 if g == YOUNG else 1 for g in groups]
    lr = LogisticRegression()
    lr.fit(x_r, y)
    p_y = [1 if x[0] < 0.5 else 0 for x in lr.predict_proba(x_r)]
    error = np.sum(np.abs(np.subtract(y, p_y)))
    ax.set_title('{} error: {}'.format(title, error))
    return error


class Plot(Enum):
    SCATTER = 1
    HIST = 2
    MA = 3


def mean_regions(df, title, ax, plot_type):
    """Plots for mean values over OD and YD"""
    ods = [c for c in df.columns if is_od(c)]
    yds = [c for c in df.columns if is_yd(c)]

    signal = pd.DataFrame()
    signal["ODS"] = df[ods].mean(axis=1)
    signal["YDS"] = df[yds].mean(axis=1)

    if plot_type == Plot.MA:
        signal["M"] = signal["ODS"] / signal["YDS"]
        # Fix division by zero to 0.0
        signal.loc[~np.isfinite(signal["M"]), "M"] = 0.0
        signal["A"] = 0.5 * (signal["ODS"] + signal["YDS"])
        ax.scatter(signal["A"], signal["M"], alpha=.3, s=1)
        ax.set_xlabel("A")
        ax.set_ylabel("M")

        xmin = np.min(ax.get_xlim())
        xmax = np.max(ax.get_xlim())
        ax.plot([xmin, xmax], [0, 0], c="red", alpha=0.75, lw=1, ls='dotted')
        ax.set_xlim([xmin, xmax])

    elif plot_type == Plot.HIST:
        signal["ODS"] = signal["ODS"]
        signal["YDS"] = signal["YDS"]

        ax.hist(signal["ODS"], color=OLD.color, bins=100, alpha=0.3, label="ODS")
        ax.hist(signal["YDS"], color=YOUNG.color, bins=100, alpha=0.3, label="YDS")
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
        # Line x = y
        ax.plot(lims, lims, 'r-', alpha=0.75, lw=1, ls='dotted')
        ax.set_aspect('equal')
        ax.set_xlim(lims)
        ax.set_ylim(lims)

    ax.set_title(title)


def mean_boxplots(df, title, ax):
    """Plot mean values for individual donors"""
    signal = df.mean(axis=1).to_frame("value")
    signal.index = [age(n) for n in signal.index]
    # Setup age
    signal["age"] = "ODS"
    signal.loc[[bool(is_yd(str(x))) for x in signal.index], "age"] = "YDS"

    age_labels = list(reversed(sorted(set(signal['age']))))
    sns.boxplot(x="age", y="value", data=signal, palette="Set3",
                linewidth=1.0, order=age_labels, ax=ax)
    sns.swarmplot(x="age", y="value", data=signal, color=".25", order=age_labels, ax=ax)

    for i, age_label in enumerate(age_labels):
        age_data = signal[signal['age'] == age_label]
        for j, label in enumerate(age_data.index):
            ax.annotate(label, xy=(i, age_data.iloc[j, :]['value']),
                        xytext=(5, 0),
                        color=OLD.color if age_label == "YDS" else YOUNG.color,
                        textcoords='offset points')

    ax.set_title(title)


def visualize(f, signal_type):
    try:
        print('Visualizing', signal_type, f)
        df = pd.read_csv(f, sep='\t')
        od_inputs = [c for c in df.columns.values if is_od_input(c)]
        yd_inputs = [c for c in df.columns.values if is_yd_input(c)]
        if od_inputs and yd_inputs:
            signal = df.drop(['chr', 'start', 'end', od_inputs[0], yd_inputs[0]], axis=1)
        else:
            signal = df.drop(['chr', 'start', 'end'], axis=1)

        plt.figure(figsize=(30, 6))
        fit_error = signal_pca(signal.T, title=signal_type, ax=plt.subplot(1, 5, 1))
        mean_regions(df, title=signal_type, ax=plt.subplot(1, 5, 2),
                     plot_type=Plot.SCATTER)
        mean_regions(df, title='MA {}'.format(signal_type), ax=plt.subplot(1, 5, 3),
                     plot_type=Plot.MA)
        mean_regions(df, title='Log {}'.format(signal_type), ax=plt.subplot(1, 5, 4),
                     plot_type=Plot.HIST)
        mean_boxplots(signal.T, title=signal_type, ax=plt.subplot(1, 5, 5))
        plt.savefig(re.sub('.tsv', '.png', f))
        plt.close()

        # Save pca fit errors to file
        pd.DataFrame(data=[[fit_error]]).to_csv(re.sub('.tsv', '_pca_fit_error.csv', f),
                                                index=None, header=False)

    except FileNotFoundError as e:
        print(e)


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
