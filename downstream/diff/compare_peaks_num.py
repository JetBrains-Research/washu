import argparse
import os
import re
import subprocess
from pathlib import Path
import numpy as np
from scipy.stats import fisher_exact

import pandas as pd

import downstream.loci_of_interest as loi
from bed.bedtrace import consensus
from downstream.consensus import save_cons_to_file

help_data = """
Usage:
    compare_peaks_num.py [input folder] [outliers path] [output folder] 

Scripts search for difference based on number of track with peaks 
"""


def p_adjust_bh(p):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]


def _cli():
    parser = argparse.ArgumentParser(description=help_data,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("peaks", help="Peaks folder")
    parser.add_argument("outliers_path", help="Outlier file path")
    parser.add_argument("output", help="Output folder")

    args = parser.parse_args()
    peaks_path = Path(args.peaks)
    output_path = Path(args.output)
    outliers_path = args.outliers_path
    outliers_df = pd.read_csv(outliers_path, delimiter="\t", skiprows=1, index_col="donor")

    for hist_mod in ["H3K27ac", "H3K27me3", "H3K36me3", "H3K4me3", "H3K4me1"]:
        tool = "zinbra"

        print("Processing {}".format(hist_mod))
        tool_path = peaks_path / hist_mod / tool
        tracks_paths = sorted({tool_path / file for file in os.listdir(str(tool_path)) if
                               re.match('.*(?:_peaks.bed)$', file)},
                              key=loi.donor_order_id)
        tracks_names = list({str(tracks_path) for tracks_path in tracks_paths})
        od_paths_map = {re.findall('OD\\d+', track_name)[0]:
                            track_name for track_name in tracks_names if re.match('.*OD\\d+.*',
                                                                                  track_name)}
        yd_paths_map = {re.findall('YD\\d+', track_name)[0]:
                            track_name for track_name in tracks_names if re.match('.*YD\\d+.*',
                                                                                  track_name)}

        for donor in outliers_df.loc[:, hist_mod][outliers_df.loc[:, hist_mod] == 1].keys():
            if donor in od_paths_map.keys():
                del od_paths_map[donor]
            if donor in yd_paths_map.keys():
                del yd_paths_map[donor]

        weak_cons = consensus(list(od_paths_map.values()) + list(yd_paths_map.values()), 0, 30)
        consensus_path = output_path / "{}_consensus.bed".format(hist_mod)
        save_cons_to_file(weak_cons, consensus_path)

        find_diff(list(od_paths_map.values()), list(yd_paths_map.values()),
                  output_path, hist_mod + "_" + tool, consensus_path)


def find_diff(od_values, yd_values, output_path, prefix, consensus_path):
    merged_path = "{}/{}_result.bed".format(output_path, prefix)

    y_peaks = count_peaks(od_values, output_path, "OD", consensus_path)
    o_peaks = count_peaks(yd_values, output_path, "YD", consensus_path)

    cmd = "cat {} | ".format(" ".join(y_peaks + o_peaks))
    cmd += "sort -k 1,1 -k2,2n | bedtools merge -c 4 -o collapse >{}".format(merged_path)
    subprocess.check_call(["bash", "-c", cmd])

    for file in y_peaks + o_peaks:
        os.remove(file)

    on = len(od_values)
    yn = len(yd_values)

    lines = []
    p_vals = []
    change = []

    with open(merged_path) as f:
        for l in f:
            parts = l.split()
            op = 0
            yp = 0

            for v in parts[3].split(","):
                if v == "YD":
                    yp += 1
                if v == "OD":
                    op += 1

            oddsratio, pvalue = fisher_exact([[yp, yn - yp], [op, on - op]])
            lines.append("\t".join(parts[0:3]) + "\t{},{}".format(yp, op))
            p_vals.append(pvalue)
            change.append(op / on - yp / yn)

    os.remove(merged_path)

    p_adjs = p_adjust_bh(p_vals)

    sorted_values = sorted(zip(p_vals, p_adjs, change, lines))

    with open("{}/{}_all.bed".format(output_path, prefix), "w") as f:
        for p_val, p_adj, change, l in sorted_values:
            f.write("{}\t{}\t{}\n".format(l.strip(), p_val, p_adj))

    counter = 1

    with open("{}/{}_up.bed".format(output_path, prefix), "w") as f:
        for p_val, p_adj, change, l in sorted_values:
            if change < 0:
                continue
            f.write("{}\t{}\t{}\n".format(l.strip(), p_val, p_adj))
            counter += 1

            if counter > 500:
                break

    counter = 1

    with open("{}/{}_down.bed".format(output_path, prefix), "w") as f:
        for p_val, p_adj, change, l in sorted_values:
            if change > 0:
                continue
            f.write("{}\t{}\t{}\n".format(l.strip(), p_val, p_adj))
            counter += 1

            if counter > 500:
                break


def count_peaks(files, output_path, prefix, consensus_path):
    intersected_files = []
    for file in files:
        intersection_file = output_path / os.path.basename(file)
        cmd = "bedtools intersect -wa -u -a {} -b {} | sed \'s/$/\t{}/\' >{}".format(
            consensus_path, file, prefix, intersection_file)
        subprocess.check_call(["bash", "-c", cmd])
        intersected_files.append(str(intersection_file))

    return intersected_files


if __name__ == "__main__":
    _cli()
