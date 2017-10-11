import os
import sys
import pybedtools
import multiprocessing
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import re
import functools
from matplotlib.backends.backend_pdf import PdfPages
from pybedtools.helpers import BEDToolsError

from bed.bedtrace import intersect, Bed, metapeaks

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
    plt.xticks(ind, names, rotation=70)
    plt.legend((p1[0], p2[0], p3[0]), ('Common', 'Group', 'Individual'))
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
    tracks_names = [x.split('/')[-1] for x in sorted(list(set(tracks_paths)))]
    callers_for_heatmaps = [item.split('_')[0] + "_" + item.split('_')[1] for item in tracks_names]
    # callers_for_heatmaps = tracks_names
    help_dict = {tracks_names[n]: n for n in range(len(tracks_names))}
    heatmap = np.zeros((len(tracks_names), len(tracks_names)))
    sample = []

    for file_path in tracks_paths:
        distances = list(pool.map(functools.partial(calc_jaccard_distance, bed_b=file_path,
                                                    size=len(tracks_names) * len(tracks_names)), tracks_paths))
        short_distances = [(y, x) for (y, x) in sorted(zip(distances, tracks_names), reverse=True)]
        for db_entry in short_distances:
            heatmap[help_dict[file_path.split('/')[-1]], help_dict[db_entry[1]]] = db_entry[0]
            sample.append(db_entry[0])

    sample = sorted(sample)
    plt.figure()
    plt.plot(sample, np.arange(len(sample)) / len(sample), label="")
    plt.xlabel('Jaccard similarity')
    plt.ylabel('CDF')
    plt.legend(loc=4)
    plt.tight_layout()
    pp.savefig()

    plt.figure()
    plt.imshow(heatmap, aspect='auto', cmap="viridis", interpolation="nearest")
    plt.xticks(np.arange(0, len(tracks_names), 1), callers_for_heatmaps, rotation='vertical')
    plt.yticks(np.arange(0, len(tracks_names), 1), callers_for_heatmaps)
    plt.colorbar(orientation='vertical')
    plt.tight_layout()
    plt.grid()
    pp.savefig()


def calc_jaccard_distance(bed_a, bed_b, size):
    with lock:
        counter.value += 1
        if counter.value % 100 == 0:
            print("{0} / {1} Processing file {2}.".format(counter.value, size, bed_a))
            pybedtools.cleanup()

    a = pybedtools.bedtool.BedTool(bed_a).sort()
    b = pybedtools.bedtool.BedTool(bed_b).sort()
    try:
        jac_dict = pybedtools.bedtool.BedTool.jaccard(a, b)
    except BEDToolsError:
        print('Failed to compare {} and {}.'.format(bed_a, bed_b))
        jac_dict = {'jaccard': -1}

    return jac_dict['jaccard']


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
    plt.figure()
    bins = np.logspace(2.0, 4.0, 80)
    for i, track_path in enumerate(tracks_paths):
        lenghts = []

        for line in read_all_lines(track_path):
            parts = line.split("\t")
            lenghts.append(int(parts[2]) - int(parts[1]))

        ax = plt.subplot(330 + i % 9 + 1)
        ax.hist(lenghts, bins, histtype='bar')
        ax.set_xscale('log')
        ax.set_xlabel('Peaks length')
        ax.set_ylabel('Peaks count')
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
                              re.match('.*peaks.*\.(?:broadPeak|bed|narrowPeak)$', f)})
    tracks_paths = sorted({bed_file for bed_file in bed_files_paths if re.match(".*([YO]D\d+).*", bed_file)})
    od_paths_map = {re.findall('OD\\d+', track_path)[0]: Bed(track_path) for track_path in tracks_paths
                    if re.match('.*OD\\d+.*\.(?:broadPeak|bed|narrowPeak)$', track_path)}
    yd_paths_map = {re.findall('YD\\d+', track_path)[0]: Bed(track_path) for track_path in tracks_paths
                    if re.match('.*YD\\d+.*\.(?:broadPeak|bed|narrowPeak)$', track_path)}
    rip_files = sorted([folder_path + '/' + file for file in os.listdir(folder_path) if file.endswith("_rip.csv")])

    pp = PdfPages(args[2])

    plot_peaks_intersect(od_paths_map, yd_paths_map, pp)
    plot_consensus(tracks_paths, pp)
    plot_jaccard_heatmap(bed_files_paths, pp)
    plot_frip_boxplot(rip_files, pp)
    plot_length_hist(tracks_paths, pp)

    pp.close()


if __name__ == "__main__":
    num_of_threads = 20
    # Initializing multiprocessor pool
    counter = multiprocessing.Value('i', 0)
    lock = multiprocessing.Lock()
    pool = multiprocessing.Pool(processes=num_of_threads)

    main()
