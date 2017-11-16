import functools
import tempfile
import operator
import multiprocessing
import numpy as np
import pandas as pd
from pathlib import Path

# Force matplotlib to not use any Xwindows backend.
import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt  # nopep8
import seaborn as sns  # nopep8
from scripts.util import age, is_od_or_yd  # nopep8
from bed.bedtrace import intersect, Bed, metapeaks, union  # nopep8
from reports.bed_metrics import save_plot  # nopep8


def venn_consensus(od_consensus_bed, yd_consensus_bed, scale, save_to=None):
    """
    Plots venn diagram for consensus of selected scale:

    :param od_consensus_bed: OD bed with od group consensus
    :param yd_consensus_bed: YD bed with yd group consensus
    :param scale: 1/scale is ratio of tack number needed for consensus
    :param save_to: Object for plots saving
    """
    plt.figure()
    plt.title("Required consensus: %.2f%%" % (100.0 / scale))
    metapeaks({'Young donors': yd_consensus_bed, 'Old donors': od_consensus_bed})
    save_plot(save_to)


def bar_consensus(od_paths_map, yd_paths_map, od_consensus_bed, yd_consensus_bed, yd_od_int_bed,
                  threads_num, save_to=None):
    """
    Plots venn diagram and bar plot for consensus of selected scale:

    :param od_paths_map: OD names as keys, od beds as values
    :param yd_paths_map: YD names as keys, yd beds as values
    :param od_consensus_bed: OD bed with od group consensus
    :param yd_consensus_bed: YD bed with yd group consensus
    :param yd_od_int_bed: BED with intersection of od and yd groups consensuses
    :param save_to: Object for plots saving
    :param threads_num: Threads number for parallel execution
    """
    pool = multiprocessing.Pool(processes=threads_num)
    n = len(yd_paths_map) + len(od_paths_map)
    ind = np.arange(n)

    yd_names, yd_common, yd_own, yd_opposite, yd_personal = zip(*pool.map(
        functools.partial(groups_sizes, common_bed=yd_od_int_bed,
                          own_group_bed=yd_consensus_bed,
                          opposite_group_bed=od_consensus_bed), sorted(yd_paths_map.items(),
                                                                       key=operator.itemgetter(0))))
    od_names, od_common, od_own, od_opposite, od_personal = zip(*pool.map(
        functools.partial(groups_sizes, common_bed=yd_od_int_bed,
                          own_group_bed=od_consensus_bed,
                          opposite_group_bed=yd_consensus_bed), sorted(od_paths_map.items(),
                                                                       key=operator.itemgetter(0))))

    common_peaks = yd_common + od_common
    group_own = yd_own + od_own
    group_opposite = yd_opposite + od_opposite

    plt.figure()
    width = 0.35
    p1 = plt.bar(ind, common_peaks, width, color='green')
    p2 = plt.bar(ind, group_own, width, bottom=[yd_od_int_bed.count()] * n,
                 color='blue')
    p3 = plt.bar(ind, group_opposite, width, bottom=[yd_od_int_bed.count() +
                                                     max(group_own)] * n, color='orange')
    p4 = plt.bar(ind, yd_personal + od_personal, width,
                 bottom=[yd_od_int_bed.count() + max(group_own) +
                         max(group_opposite)] * n, color='black')
    plt.ylabel('Peaks count')
    plt.xticks(ind, yd_names + od_names, rotation=90)
    plt.legend((p1[0], p2[0], p3[0], p4[0]),
               ('Common', 'Own Group', 'Opposite Group', 'Individual'))
    plt.tight_layout()
    save_plot(save_to)


def calc_consensus(od_paths_map, yd_paths_map, scale):
    """
    Calculates consensus of selected scale:

    :param od_paths_map: OD names as keys, od beds as values
    :param yd_paths_map: YD names as keys, yd beds as values
    :param scale: 1/scale is ratio of tack number needed for consensus
    :return OD consensus bed, YD consensus bed and bed with their intersect
    """
    od_union = union(*[path[1] for path in sorted(od_paths_map.items(),
                                                  key=operator.itemgetter(0))])
    yd_union = union(*[path[1] for path in sorted(yd_paths_map.items(),
                                                  key=operator.itemgetter(0))])
    od_union.compute()
    yd_union.compute()

    od_consensus = [region for region in od_union.cat().split('\n') if region.count('|') >=
                    np.ceil(float(len(od_paths_map)) / scale) - 1.0]
    yd_consensus = [region for region in yd_union.cat().split('\n') if region.count('|') >=
                    np.ceil(float(len(yd_paths_map)) / scale) - 1.0]

    od_consensus_path = tempfile.NamedTemporaryFile(mode='w', suffix='.bed', prefix='od_consensus',
                                                    delete=False)
    yd_consensus_path = tempfile.NamedTemporaryFile(mode='w', suffix='.bed', prefix='yd_consensus',
                                                    delete=False)
    od_consensus_path.write('\n'.join(od_consensus))
    yd_consensus_path.write('\n'.join(yd_consensus))
    od_consensus_path.close()
    yd_consensus_path.close()

    od_consensus_bed = Bed(od_consensus_path.name)
    yd_consensus_bed = Bed(yd_consensus_path.name)
    yd_od_int_bed = intersect(od_consensus_bed, yd_consensus_bed)
    yd_od_int_bed.compute()

    return od_consensus_bed, yd_consensus_bed, yd_od_int_bed


def groups_sizes(entry, common_bed, own_group_bed, opposite_group_bed):
    """
    Count sizes of common, own group, opposite group and personal peaks:

    :param entry: Pair of donor name and bed file
    :param common_bed: Bed file with common peaks
    :param own_group_bed: Bed file with own group peaks
    :param opposite_group_bed: Bed file opposite group peaks
    :return Donor name and number of peaks in common, own group, opposite group and personal
    """
    common = intersect(entry[1], common_bed).count()
    own_group = intersect(entry[1], own_group_bed).count() - common
    opposite_group = intersect(entry[1], opposite_group_bed).count() - common
    personal = entry[1].count() - common - own_group - opposite_group
    return entry[0], common, own_group, opposite_group, personal


def cumulative_consensus(tracks_paths, save_to=None):
    """
    Plot cumulative consensus from individual number of peaks down to common number of peaks:

    :param tracks_paths: List of tracks paths
    :param save_to: Object for plots saving
    """
    tracks_union = union(*[Bed(track_path) for track_path in tracks_paths])
    tracks_union.compute()

    counts = [0] * len(tracks_paths)
    for line in tracks_union.cat().split('\n'):
        if line != '':
            parts = line.split("\t")
            count = len(parts[3].split("|"))
            counts[count - 1] += 1
    counts.reverse()
    counts_cumulative = list(np.cumsum(counts))
    counts_cumulative.reverse()

    plt.figure()
    plt.xlabel('Number of donors')
    plt.ylabel('Number of peaks')
    plt.plot(range(1, len(tracks_paths) + 1), counts_cumulative, label="")
    plt.title("Reverse cumulative consensus peaks sum via number of donors")
    save_plot(save_to)


def frip_peaks(age_labels, df, save_to):
    """
    Plots FRiP via peaks number for passed in data frame donors:

    :param age_labels: Age labels for dots coloring
    :param df: Data frame with information about donors and their FRiP
    :param save_to: Object for plots saving
    """
    plt.figure()
    ax = plt.subplot()
    for i, age_label in enumerate(age_labels):
        age_data = df[df['age'] == age_label]
        plt.plot(age_data["peaks"], age_data["frip"], 'ro',
                 color="red" if age_label == "YDS" else "blue", label=age_label)
        for j, label in enumerate(age_data.index):
            ax.annotate(label, xy=(age_data["peaks"][j], age_data["frip"][j]), xytext=(5, 0),
                        color="red" if age_label == "YDS" else "blue", textcoords='offset points')
    plt.xlabel('Peaks')
    plt.ylabel('FRiP')
    plt.legend(loc=4)
    plt.tight_layout()
    save_plot(save_to)


def frip_boxplot(age_labels, df, save_to):
    """
    Plots FRiP boxplot for passed in data frame donors:

    :param age_labels: Age labels for dots coloring
    :param df: Data frame with information about donors and their FRiP
    :param save_to: Object for plots saving
    """
    plt.figure()
    ax = plt.subplot()
    sns.boxplot(x="age", y="frip", data=df, palette="Set3", linewidth=1.0, order=age_labels, ax=ax)
    sns.swarmplot(x="age", y="frip", data=df, color=".25", order=age_labels, ax=ax)

    for i, age_label in enumerate(age_labels):
        age_data = df[df['age'] == age_label]
        for j, label in enumerate(age_data.index):
            ax.annotate(label, xy=(i, age_data.iloc[j, :]['frip']), xytext=(5, 0),
                        color="red" if age_label == "YDS" else "blue", textcoords='offset points')
    ax.set_title("Signal FRiP")
    save_plot(save_to)


def calc_frip(rip_paths):
    """
    Calculates FRiP for old and young donor groups:

    :param rip_paths: List of absolute paths to rip files
    :return Age labels and data frame with information about donors FRiP
    """
    df = pd.DataFrame(columns=["file", "reads", "peaks", "rip"])
    for rip_path in rip_paths:
        data = pd.read_csv(rip_path)
        rip_file = rip_path.split('/')[-1]
        reads = float(data["reads"].loc[0])
        peaks = int(data["peaks"].loc[0])
        rip = float(data["rip"].loc[0])
        df.loc[df.size] = (rip_file, reads, peaks, rip)
    df["frip"] = 100.0 * df["rip"] / df["reads"]
    df.index = [age(df.loc[n]["file"]) for n in df.index]
    df["age"] = "ODS"
    df.loc[df.index.str.startswith("Y"), "age"] = "YDS"
    age_labels = list(reversed(sorted(list(set(df['age'])))))
    return age_labels, df


def length_bar_plots(tracks_paths, min_power, max_power, threads_num, save_to=None):
    """
    Plots bar plot for each track - peaks count via peaks lengths:

    :param tracks_paths: List of absolute track paths
    :param min_power: used for left border of bar plot as a power for 10
    :param max_power: used for right border of bar plot as a power for 10
    :param threads_num: Threads number for parallel execution
    :param save_to: Object for plots saving
    """
    pool = multiprocessing.Pool(processes=threads_num)
    bins = np.logspace(min_power, max_power, 80)
    ordered_paths, ordered_lengths, track_max_bar_height = zip(*pool.map(functools.partial(
        _calculate_lengths, bins=bins), tracks_paths))
    max_bar_height = max(track_max_bar_height)
    lengths = dict(zip(ordered_paths, ordered_lengths))

    plt.figure()
    for i, track_path in enumerate(tracks_paths):
        ax = plt.subplot(330 + i % 9 + 1)
        ax.hist(lengths[track_path], bins, histtype='bar')
        ax.set_xscale('log')
        ax.set_xlabel('Peaks length')
        ax.set_ylabel('Peaks count')
        ax.set_ylim([0, max_bar_height])
        ax.set_title(age(track_path) if is_od_or_yd(track_path) else Path(track_path).name)

        if i % 9 == 8:
            plt.tight_layout()
            save_plot(save_to)
            plt.figure()
    plt.tight_layout()
    save_plot(save_to)


def _calculate_lengths(track_path, bins):
    """
    Calculates track peaks lengths and max bar height within all bins:

    :param track_path: Absolute track path
    :param bins: Object with information about bins in which peaks should be divided
    :return Absolute track path, list of track peaks lengths, max bar height for the track
    """
    lengths = []
    extended_bins = np.append(bins, bins[-1] + 1)
    with open(track_path, "r") as track_file:
        for line in track_file:
            parts = line.split("\t")
            lengths.append(int(parts[2]) - int(parts[1]))
    counts = [0] * len(extended_bins)
    for bin_index in np.digitize(lengths, extended_bins):
        counts[bin_index - 1] += 1
    return track_path, lengths, max(counts[:len(bins)])