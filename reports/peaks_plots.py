import os
import sys
import multiprocessing
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import re
import functools
import tempfile

from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster import hierarchy
from scipy.spatial import distance
from bed.bedtrace import intersect, Bed, metapeaks, jaccard, union, run

__author__ = 'petr.tsurinov@jetbrains.com'
help_data = """
Usage: peaks_plots.py [peaks folder] [output pdf path] [top peaks count (optional)]

Script creates pdf report with ChIP-seq peaks statistics:
 1) median peak consensus venn diagram
 2) median peak consensus bar plot
 3) cumulative consensus plot
 4) Jaccard similarity plot
 5) Jaccard index heatmap
 6) Frip/peaks plot
 7) Frip/age boxplot
 8) Peaks count/peaks length bar plots for each donor
"""


def plot_venn_bar_consensus(od_paths_map, yd_paths_map, pp, scale):
    od_union = union(*od_paths_map.values())
    yd_union = union(*yd_paths_map.values())
    od_union.compute()
    yd_union.compute()

    od_median_consensus = [region for region in od_union.cat().split('\n') if region.count('|') >=
                           np.ceil(float(len(od_paths_map)) / scale) - 1.0]
    yd_median_consensus = [region for region in yd_union.cat().split('\n') if region.count('|') >=
                           np.ceil(float(len(yd_paths_map)) / scale) - 1.0]

    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', prefix='od_median_consensus') as od_median_consensus_path:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', prefix='yd_median_consensus') \
                as yd_median_consensus_path:
            od_median_consensus_path.write('\n'.join(od_median_consensus))
            yd_median_consensus_path.write('\n'.join(yd_median_consensus))

            yd_od_intersection = intersect(Bed(od_median_consensus_path.name), Bed(yd_median_consensus_path.name))
            yd_od_intersection.compute()

            plt.figure()
            plt.title("Required consensus: %.2f%%" % (100.0 / scale))
            metapeaks({'Young donors': Bed(yd_median_consensus_path.name),
                       'Old donors': Bed(od_median_consensus_path.name)})
            pp.savefig()

            n = len(yd_paths_map) + len(od_paths_map)
            ind = np.arange(n)
            yd_median_consensus_bed = Bed(yd_median_consensus_path.name)
            od_median_consensus_bed = Bed(od_median_consensus_path.name)

            yd_names, yd_common_peaks, yd_first_group, yd_second_group_specific, yd_sample_specific = zip(
                *pool.map(functools.partial(calculate_intersects, common_bed=yd_od_intersection,
                                            own_group_bed=yd_median_consensus_bed,
                                            opposite_group_bed=od_median_consensus_bed), yd_paths_map.items()))
            od_names, od_common_peaks, od_first_group, od_second_group_specific, od_sample_specific = zip(
                *pool.map(functools.partial(calculate_intersects, common_bed=yd_od_intersection,
                                            own_group_bed=od_median_consensus_bed,
                                            opposite_group_bed=yd_median_consensus_bed), od_paths_map.items()))

            common_peaks = yd_common_peaks + od_common_peaks
            first_group_specific = yd_first_group + od_first_group
            second_group_specific = yd_second_group_specific + od_second_group_specific

            plt.figure()
            width = 0.35
            p1 = plt.bar(ind, common_peaks, width, color='green')
            p2 = plt.bar(ind, first_group_specific, width, bottom=[yd_od_intersection.count()] * n, color='blue')
            p3 = plt.bar(ind, second_group_specific, width, bottom=[yd_od_intersection.count() +
                                                                    max(first_group_specific)] * n, color='orange')
            p4 = plt.bar(ind, yd_sample_specific + od_sample_specific, width, bottom=[yd_od_intersection.count() +
                                                                                      max(first_group_specific) +
                                                                                      max(second_group_specific)] * n,
                         color='black')
            plt.ylabel('Peaks count')
            plt.xticks(ind, yd_names + od_names, rotation=90)
            plt.legend((p1[0], p2[0], p3[0], p4[0]), ('Common', 'Own Group', 'Opposite Group', 'Individual'))
            plt.tight_layout()
            pp.savefig()


def calculate_intersects(entry, common_bed, own_group_bed, opposite_group_bed):
    common = intersect(entry[1], common_bed).count()
    own_group = intersect(entry[1], own_group_bed).count() - common
    opposite_group = intersect(entry[1], opposite_group_bed).count() - common
    return entry[0], common, own_group, opposite_group, entry[1].count() - common - own_group - opposite_group


def plot_cumulative_consensus(tracks_paths, pp):
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
    plt.plot(range(len(tracks_paths)), counts_cumulative, label="")
    plt.title("Reverse cumulative consensus peaks sum via number of donors")
    pp.savefig()


def plot_jaccard_heatmap(tracks_paths, pp, suffix="_caller"):
    callers_for_heatmaps = [re.match(".*([YO]D\d+|consensus).*", item).group(1) + suffix for item in tracks_paths]
    help_dict = {tracks_paths[n]: n for n in range(len(tracks_paths))}
    heatmap = np.zeros((len(tracks_paths), len(tracks_paths)))

    sample = []
    for file_path in tracks_paths:
        distances = list(pool.map(functools.partial(calc_jaccard_distance, bed_b=file_path,
                                                    size=len(tracks_paths) * len(tracks_paths)), tracks_paths))
        short_distances = [(y, x) for (y, x) in sorted(zip(distances, tracks_paths), reverse=True)]
        for db_entry in short_distances:
            heatmap[help_dict[file_path], help_dict[db_entry[1]]] = db_entry[0]
            sample.append(db_entry[0])

    sample = sorted(sample)
    plt.figure()
    plt.plot(sample, np.arange(len(sample)) / len(sample), label="")
    plt.xlabel('Jaccard similarity')
    plt.ylabel('CDF')
    plt.legend(loc=4)
    plt.tight_layout()
    pp.savefig()

    dissimilarity = distance.squareform(1 - heatmap)
    linkage = hierarchy.linkage(dissimilarity, method="average")
    clusters = hierarchy.dendrogram(linkage, no_plot=True)['leaves']
    heatmap_data_frame = pd.DataFrame(data=heatmap, index=callers_for_heatmaps, columns=callers_for_heatmaps)

    figure, axes = plt.subplots()
    sns.set(font_scale=0.75)
    cg = sns.clustermap(heatmap_data_frame, cmap='rainbow', row_linkage=linkage, col_linkage=linkage,
                        row_cluster=clusters, col_cluster=clusters, vmin=0.0, vmax=1.0)
    plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
    figure.tight_layout()
    pp.savefig()


def calc_jaccard_distance(bed_a, bed_b, size):
    with lock:
        counter.value += 1
        if counter.value % 100 == 0:
            print("{0} / {1} Processing file {2}.".format(counter.value, size, bed_a))
    return 1.0 if bed_a == bed_b else jaccard(bed_a, bed_b)


def plot_frip_boxplot(rip_files, pp):
    df = pd.DataFrame(columns=["file", "reads", "peaks", "rip"])
    for rip_file in rip_files:
        data = pd.read_csv(rip_file)
        rip_file = rip_file.split('/')[-1]
        reads = float(data["reads"].loc[0])
        peaks = int(data["peaks"].loc[0])
        rip = float(data["rip"].loc[0])
        df.loc[df.size] = (rip_file, reads, peaks, rip)
    df["frip"] = 100.0 * df["rip"] / df["reads"]
    df.index = [re.search('[yo]d\\d+', df.loc[n]["file"], flags=re.IGNORECASE).group(0) for n in df.index]
    df["age"] = "ODS"
    df.loc[df.index.str.startswith("Y"), "age"] = "YDS"
    age_labels = list(reversed(sorted(list(set(df['age'])))))

    plt.figure()
    ax = plt.subplot()
    for i, age_label in enumerate(age_labels):
        age_data = df[df['age'] == age_label]
        plt.plot(age_data["peaks"], age_data["frip"], 'ro', color="red" if age_label == "YDS" else "blue",
                 label=age_label)
        for j, label in enumerate(age_data.index):
            ax.annotate(label, xy=(age_data["peaks"][j], age_data["frip"][j]), xytext=(5, 0),
                        color="red" if age_label == "YDS" else "blue", textcoords='offset points')
    plt.xlabel('Peaks')
    plt.ylabel('FRiP')
    plt.legend(loc=4)
    plt.tight_layout()
    pp.savefig()

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
    pp.savefig()


def plot_length_bar(tracks_paths, pp, min_length_power, max_length_power):
    bins = np.logspace(min_length_power, max_length_power, 80)
    track_path, lengths, track_max_bar_height = zip(*pool.map(functools.partial(calculate_lengths, bins=bins),
                                                              tracks_paths))
    max_length = max(track_max_bar_height)
    lengths = dict(zip(track_path, lengths))

    plt.figure()
    for i, track_path in enumerate(tracks_paths):
        ax = plt.subplot(330 + i % 9 + 1)
        ax.hist(lengths[track_path], bins, histtype='bar')
        ax.set_xscale('log')
        ax.set_xlabel('Peaks length')
        ax.set_ylabel('Peaks count')
        ax.set_ylim([0, max_length])
        ax.set_title(re.match(".*([YO]D\d+).*", track_path).group(1))

        if i % 9 == 8:
            plt.tight_layout()
            pp.savefig()
            plt.figure()
    plt.tight_layout()
    pp.savefig()


def calculate_lengths(track_path, bins):
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


def main():
    args = sys.argv

    if len(args) < 2:
        print(help_data)
        sys.exit(1)

    folder_path = args[1]
    bed_files_paths = sorted({folder_path + '/' + f for f in os.listdir(folder_path) if
                              re.match('.*\.(?:broadPeak|bed|narrowPeak)$', f)})
    filtered_bed_files_paths = []
    if len(args) == 4:
        for bed_files_path in bed_files_paths:
            with open(tempfile.tempdir + "/" + bed_files_path.split("/")[-1].split(".")[0] + "_" + args[3] + ".bed",
                      'w') as tmpfile:
                run((["sort", "-k9nr", bed_files_path], ["head", "-n", args[3]]), stdout=tmpfile)
                filtered_bed_files_paths.append(tmpfile.name)
    else:
        filtered_bed_files_paths = bed_files_paths

    tracks_paths = sorted({bed_file for bed_file in filtered_bed_files_paths if re.match(".*([YO]D\d+).*", bed_file)})
    od_paths_map = {re.findall('OD\\d+', track_path)[0]: Bed(track_path) for track_path in tracks_paths
                    if re.match('.*OD\\d+.*\.(?:broadPeak|bed|narrowPeak)$', track_path)}
    yd_paths_map = {re.findall('YD\\d+', track_path)[0]: Bed(track_path) for track_path in tracks_paths
                    if re.match('.*YD\\d+.*\.(?:broadPeak|bed|narrowPeak)$', track_path)}
    rip_files = sorted([folder_path + '/' + file for file in os.listdir(folder_path) if file.endswith("_rip.csv")])

    pp = PdfPages(args[2])

    # Code for different consensuses investigation
    # for scale in [1, 1.1667, 1.25, 1.5, 2, 3, 5, 7]:
    #     plot_venn_bar_consensus(od_paths_map, yd_paths_map, pp, scale)

    print("Calculating median consensus")
    plot_venn_bar_consensus(od_paths_map, yd_paths_map, pp, 2.0)
    print("Calculating cumulative consensus")
    plot_cumulative_consensus(tracks_paths, pp)
    print("Calculating jaccard indexes")
    plot_jaccard_heatmap(filtered_bed_files_paths, pp, "_zinbra")
    print("Calculating frip vs age")
    plot_frip_boxplot(rip_files, pp)
    print("Calculating peaks count vs length")
    plot_length_bar(tracks_paths, pp, 2.0, 4.0)

    pp.close()


if __name__ == "__main__":
    num_of_threads = 8
    # Initializing multiprocessor pool
    counter = multiprocessing.Value('i', 0)
    lock = multiprocessing.Lock()
    pool = multiprocessing.Pool(processes=num_of_threads)

    main()
