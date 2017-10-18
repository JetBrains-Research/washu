import os
import sys
import multiprocessing
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import re
import functools
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster import hierarchy
from scipy.spatial import distance

from bed.bedtrace import intersect, Bed, metapeaks, jaccard

help_data = """
"""


def plot_peaks_intersect(od_paths_map, yd_paths_map, pp):
    yd_intersection = intersect(*yd_paths_map.values())
    od_intersection = intersect(*od_paths_map.values())
    yd_od_intersection = intersect(yd_intersection, od_intersection)
    metapeaks({'Young donors': yd_intersection, 'Old donors': od_intersection})
    pp.savefig()

    n = len(yd_paths_map) + len(od_paths_map)
    ind = np.arange(n)

    common_peaks = [yd_od_intersection.count()] * n
    group_specific = [yd_intersection.count() - yd_od_intersection.count()] * len(yd_paths_map) + \
                     [od_intersection.count() - yd_od_intersection.count()] * len(od_paths_map)
    sample_specific = []
    names = []
    for k, v in yd_paths_map.items():
        sample_specific.append(v.count() - yd_intersection.count())
        names.append(k)
    for k, v in od_paths_map.items():
        sample_specific.append(v.count() - od_intersection.count())
        names.append(k)

    plt.figure()
    width = 0.35
    p1 = plt.bar(ind, common_peaks, width, color='green')
    p2 = plt.bar(ind, group_specific, width, bottom=common_peaks, color='blue')
    p3 = plt.bar(ind, sample_specific, width, bottom=np.sum([common_peaks, group_specific], axis=0), color='red')
    plt.ylabel('Peaks count')
    plt.xticks(ind, names, rotation=90)
    plt.legend((p1[0], p2[0], p3[0]), ('Common', 'Group', 'Individual'))
    plt.tight_layout()
    pp.savefig()


def plot_consensus(tracks_paths, pp):
    union_sh = os.path.realpath(os.path.join(os.path.dirname(__file__), '../bed/union.sh'))
    command = "bash {} {} >{}".format(union_sh, " ".join(tracks_paths), os.path.join("report", "counts.bed"))
    os.system(command)
    counts = [0] * len(tracks_paths)
    for line in read_all_lines(os.path.join("report", "counts.bed")):
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


def plot_jaccard_heatmap(tracks_paths, pp):
    callers_for_heatmaps = [re.match(".*([YO]D\d+|consensus).*", item).group(1) + "_zinbra" for item in tracks_paths]
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

    figure, axes = plt.subplots()
    sns.set(font_scale=0.75)

    dissimilarity = distance.squareform(1 - heatmap)
    linkage = hierarchy.linkage(dissimilarity, method="average")
    clusters = hierarchy.dendrogram(linkage, no_plot=True)['leaves']

    heatmap_data_frame = pd.DataFrame(data=heatmap,  # values
                                      index=callers_for_heatmaps,  # 1st column as index
                                      columns=callers_for_heatmaps)  # 1st row as the column names
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
            ax.annotate(label, xy=(age_data["peaks"][j], age_data["frip"][j]),
                        xytext=(5, 0),
                        color="red" if age_label == "YDS" else "blue",
                        textcoords='offset points')
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
            ax.annotate(label, xy=(i, age_data.iloc[j, :]['frip']),
                        xytext=(5, 0),
                        color="red" if age_label == "YDS" else "blue",
                        textcoords='offset points')

    ax.set_title("Signal FRiP")
    pp.savefig()


def plot_length_hist(tracks_paths, pp):
    bins = np.logspace(2.0, 4.0, 80)
    lenghts = {}
    max_length = 0
    for track_path in tracks_paths:
        lenghts[track_path] = []
        for line in read_all_lines(track_path):
            parts = line.split("\t")
            lenghts[track_path].append(int(parts[2]) - int(parts[1]))
        max_length = max(max_length, max(plt.hist(lenghts[track_path], bins, histtype='bar')[0]))

    plt.figure()
    for i, track_path in enumerate(tracks_paths):
        ax = plt.subplot(330 + i % 9 + 1)
        ax.hist(lenghts[track_path], bins, histtype='bar')
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


def read_all_lines(peak_file):
    with open(peak_file, "r") as f:
        return f.readlines()


def main():
    args = sys.argv

    if len(args) < 2:
        print(help_data)
        sys.exit(1)

    folder_path = args[1]
    bed_files_paths = sorted({folder_path + '/' + f for f in os.listdir(folder_path) if
                              re.match('.*\.(?:broadPeak|bed|narrowPeak)$', f)})
    filtered_bed_files_paths = bed_files_paths
    # filtered_bed_files_paths = []
    # for bed_files_path in bed_files_paths:
    #      with open(tempfile.tempdir + "/" + bed_files_path.split("/")[-1].split(".")[0] + "_" + args[3] + ".bed",
    #                'w') as tmpfile:
    #         run((["sort", "-k9nr", bed_files_path], ["head", "-n", args[3]]), stdout=tmpfile)
    #         filtered_bed_files_paths.append(tmpfile.name)

    tracks_paths = sorted({bed_file for bed_file in filtered_bed_files_paths if re.match(".*([YO]D\d+).*", bed_file)})
    od_paths_map = {re.findall('OD\\d+', track_path)[0]: Bed(track_path) for track_path in tracks_paths
                    if re.match('.*OD\\d+.*\.(?:broadPeak|bed|narrowPeak)$', track_path)}
    yd_paths_map = {re.findall('YD\\d+', track_path)[0]: Bed(track_path) for track_path in tracks_paths
                    if re.match('.*YD\\d+.*\.(?:broadPeak|bed|narrowPeak)$', track_path)}
    rip_files = sorted([folder_path + '/' + file for file in os.listdir(folder_path) if file.endswith("_rip.csv")])

    pp = PdfPages(args[2])

    plot_peaks_intersect(od_paths_map, yd_paths_map, pp)
    plot_consensus(tracks_paths, pp)
    plot_jaccard_heatmap(filtered_bed_files_paths, pp)
    plot_frip_boxplot(rip_files, pp)
    plot_length_hist(tracks_paths, pp)

    pp.close()


if __name__ == "__main__":
    num_of_threads = 8
    # Initializing multiprocessor pool
    counter = multiprocessing.Value('i', 0)
    lock = multiprocessing.Lock()
    pool = multiprocessing.Pool(processes=num_of_threads)

    main()
