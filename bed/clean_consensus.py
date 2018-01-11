import os
import re
import sys
import argparse
import pandas as pd
from pathlib import Path
import reports.loci_of_interest as loi

__author__ = 'petr.tsurinov@jetbrains.com'
help_data = """
Usage:
"""
outliers_path = "/mnt/stripe/bio/experiments/aging/Y20O20.outliers.csv"
outliers_df = pd.read_csv(outliers_path, delimiter="\t", skiprows=1, index_col="donor")


def _cli():
    parser = argparse.ArgumentParser(description=help_data,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("peaks", help="Peaks folder")
    parser.add_argument("tools", nargs='*', help="Tools folders")

    args = parser.parse_args()
    folder = Path(args.peaks)
    tools = args.tools

    for hist_mod in ["H3K27ac", "H3K27me3", "H3K36me3", "H3K4me3", "H3K4me1"]:
        for tool in tools:
            print("Processing {} {}".format(hist_mod, tool))
            tool_path = folder / hist_mod / tool
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

            for donor in outliers_df.loc[:, hist_mod].index:
                if outliers_df.loc[:, hist_mod][donor] == 1:
                    if donor in od_paths_map.keys():
                        del od_paths_map[donor]
                    if donor in yd_paths_map.keys():
                        del yd_paths_map[donor]

            median_cons = consensus(list(od_paths_map.values()) + list(yd_paths_map.values()), 0,
                                    50)
            save_cons_to_file(median_cons, tool_path / (hist_mod + "_" + tool +
                                                        "_median_consensus.bed"))
            od_median_cons = consensus(od_paths_map.values(), 0, 50)
            save_cons_to_file(od_median_cons, tool_path / (hist_mod + "_" + tool +
                                                           "_ODS_median_consensus.bed"))
            yd_median_cons = consensus(yd_paths_map.values(), 0, 50)
            save_cons_to_file(yd_median_cons, tool_path / (hist_mod + "_" + tool +
                                                           "_YDS_median_consensus.bed"))

            weak_cons = consensus(list(od_paths_map.values()) + list(yd_paths_map.values()), 2, 0)
            save_cons_to_file(weak_cons, tool_path / (hist_mod + "_" + tool +
                                                      "_weak_consensus.bed"))
            od_weak_cons = consensus(od_paths_map.values(), 2, 0)
            save_cons_to_file(od_weak_cons, tool_path / (hist_mod + "_" + tool +
                                                         "_ODS_weak_consensus.bed"))
            yd_weak_cons = consensus(yd_paths_map.values(), 2, 0)
            save_cons_to_file(yd_weak_cons, tool_path / (hist_mod + "_" + tool +
                                                         "_YDS_weak_consensus.bed"))


def save_cons_to_file(cons, path):
    file = open(path, 'wb')
    file.write(cons)
    file.close()


if __name__ == "__main__":
    parent_dir = os.path.dirname(os.path.realpath(__file__))
    project_root = os.path.abspath(os.path.join(parent_dir) + "/..")
    sys.path.insert(0, project_root)

    from bed.bedtrace import consensus  # nopep8

    _cli()
