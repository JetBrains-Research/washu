import argparse
import os
import re
from pathlib import Path

import pandas as pd

from bed.bedtrace import consensus
import downstream.loci_of_interest as loi

__author__ = 'petr.tsurinov@jetbrains.com'
help_data = """
Usage:
    consensus.py [input folder] [outliers path] [output folder] [tools list]

Script creates peaks consensuses files:
 1) Median consensus
 2) Old donor median consensus
 3) Young donor median consensus
 4) Weak consensus
 5) Old donor weak consensus
 6) Young donor weak consensus
"""


def _cli():
    parser = argparse.ArgumentParser(description=help_data,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("peaks", help="Peaks folder")
    parser.add_argument("outliers_path", help="Outlier file path")
    parser.add_argument("output", help="Output folder")
    parser.add_argument("tools", nargs='*', help="Tools folders")

    args = parser.parse_args()
    peaks_path = Path(args.peaks)
    output_path = Path(args.output)
    tools = args.tools
    outliers_path = args.outliers_path
    outliers_df = pd.read_csv(outliers_path, delimiter="\t", skiprows=1, index_col="donor")

    for hist_mod in ["H3K27ac", "H3K27me3", "H3K36me3", "H3K4me3", "H3K4me1"]:
        for tool in tools:
            print("Processing {} {}".format(hist_mod, tool))
            tool_path = peaks_path / hist_mod / tool
            tracks_paths = sorted({tool_path / file for file in os.listdir(str(tool_path)) if
                                   re.match('.*(?:_peaks.bed|-island.bed|Peak)$', file)},
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

            construct_consensuses(list(od_paths_map.values()), list(yd_paths_map.values()),
                                  output_path, hist_mod + "_" + tool, "median", 0, 50)
            construct_consensuses(list(od_paths_map.values()), list(yd_paths_map.values()),
                                  output_path, hist_mod + "_" + tool, "weak", 2, 0)


def construct_consensuses(od_values, yd_values, output_path, prefix, suffix, c, p):
    weak_cons = consensus(od_values + yd_values, c, p)
    save_cons_to_file(weak_cons, output_path / "{}_{}_consensus.bed".format(prefix, suffix))
    od_weak_cons = consensus(od_values, c, p)
    save_cons_to_file(od_weak_cons, output_path / "{}_ODS_{}_consensus.bed".format(prefix, suffix))
    yd_weak_cons = consensus(yd_values, c, p)
    save_cons_to_file(yd_weak_cons, output_path / "{}_YDS_{}_consensus.bed".format(prefix, suffix))


def save_cons_to_file(cons, path):
    file = open(str(path), 'wb')
    file.write(cons)
    file.close()


if __name__ == "__main__":
    _cli()
