import argparse
import datetime
import os
import re
import sys
from pathlib import Path

import pandas as pd

__author__ = 'petr.tsurinov@jetbrains.com'
help_data = """
Usage:
    combined_peaks_analysis_report.py [input folder] [output folder] [first tool] [second tool]

Script creates pdf reports with ChIP-seq peaks statistics:
 1) median peak consensus bar plot
 2) Metric #1 heatmap
"""
failed_tracks_path = "/mnt/stripe/bio/experiments/aging/Y20O20.failed_tracks.csv"
failed_tracks_df = pd.read_csv(failed_tracks_path, delimiter="\t", skiprows=1, index_col="donor")


def _cli():
    parser = argparse.ArgumentParser(description=help_data,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input_folder", help="Histones folder")
    parser.add_argument("output_folder", help="Output folder for pdf")
    parser.add_argument("tool1", help="First tool")
    parser.add_argument("tool2", help="Second tool")
    parser.add_argument('-p', '--threads', help="Threads number for parallel processing",
                        type=int, default=30)

    args = parser.parse_args()
    input_folder = Path(args.input_folder)
    output_folder = Path(args.output_folder)
    tool1 = args.tool1
    tool2 = args.tool2
    threads_num = args.threads

    for hist_mod in ["H3K27ac", "H3K27me3", "H3K36me3", "H3K4me3", "H3K4me1"]:
        peaks_folder1 = input_folder / hist_mod / tool1
        peaks_folder2 = input_folder / hist_mod / tool2
        today = datetime.datetime.today()
        pdf_path = str(output_folder / (hist_mod + "_" + tool1 + "_vs_" + tool2 + "_" +
                                        today.strftime("%d.%m.%Y") + ".pdf"))
        tracks_paths1 = sorted({peaks_folder1 / file for file in os.listdir(str(peaks_folder1)) if
                               re.match('.*(?:_peaks.bed|-island.bed|Peak)$', file)},
                               key=loi.donor_order_id)
        tracks_paths2 = sorted({peaks_folder2 / file for file in os.listdir(str(peaks_folder2)) if
                               re.match('.*(?:_peaks.bed|-island.bed|Peak)$', file)},
                               key=loi.donor_order_id)
        if len(tracks_paths1) > 0 and len(tracks_paths2) > 0:
            tracks_paths = tracks_paths1 + tracks_paths2

            tracks_names = list({str(tracks_path) for tracks_path in tracks_paths})
            od_paths_map = {re.findall('OD\\d+', track_name)[0] + detect_tool(track_name):
                            track_name for track_name in tracks_names if re.match('.*OD\\d+.*',
                                                                                  track_name)}
            yd_paths_map = {re.findall('YD\\d+', track_name)[0] + detect_tool(track_name):
                            track_name for track_name in tracks_names if re.match('.*YD\\d+.*',
                                                                                  track_name)}
            df = bed_metric_table(tracks_paths, tracks_paths, threads=threads_num)
            for donor in failed_tracks_df.loc[:, hist_mod].index:
                if failed_tracks_df.loc[:, hist_mod][donor] == 1:
                    _remove_donor_from_map(donor, od_paths_map)
                    _remove_donor_from_map(donor, yd_paths_map)
                    for df_index in df.index:
                        if (donor + "_") in df_index or (donor + ".") in df_index:
                            del df[df_index]
                            df = df.drop(df_index)

            with PdfPages(pdf_path) as pdf:
                print("Calculating median consensus")
                od_consensus_bed, yd_consensus_bed, yd_od_int_bed = \
                    calc_consensus_file(list(od_paths_map.values()), list(yd_paths_map.values()),
                                        percent=50)
                bar_consensus(od_paths_map, yd_paths_map, od_consensus_bed, yd_consensus_bed,
                              yd_od_int_bed, threads_num, pdf, (10, 10), 10)

                print("Calculating metric #1 indexes")
                sns.set(font_scale=0.75)
                g = plot_metric_heatmap("Intersection metric", df, figsize=(14, 14), save_to=pdf,
                                        row_label_converter=label_converter_donor_and_tool,
                                        col_label_converter=label_converter_donor_and_tool,
                                        show_or_save_plot=False)
                plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
                save_plot(pdf)

                desc = pdf.infodict()
                desc['Title'] = 'Report: Combined peaks plots for different callers'
                desc['Author'] = 'JetBrains Research BioLabs'
                desc['Subject'] = 'peaks'
                desc['CreationDate'] = datetime.datetime.today()
                desc['ModDate'] = datetime.datetime.today()


def _remove_donor_from_map(donor, paths_map):
    for path_key in list(paths_map.keys()):
        if (donor + "_") in path_key:
            del paths_map[path_key]


if __name__ == "__main__":
    # Force matplotlib to not use any Xwindows backend.
    import matplotlib

    matplotlib.use('Agg')

    parent_dir = os.path.dirname(os.path.realpath(__file__))
    project_root = os.path.abspath(os.path.join(parent_dir) + "/..")
    sys.path.insert(0, project_root)

    import seaborn as sns  # nopep8
    import downstream.loci_of_interest as loi
    import matplotlib.pyplot as plt  # nopep8
    from matplotlib.backends.backend_pdf import PdfPages  # nopep8
    from downstream.bed_metrics import bed_metric_table, plot_metric_heatmap, \
        label_converter_donor_and_tool, save_plot  # nopep8
    from downstream.peak_metrics import calc_consensus_file, bar_consensus, detect_tool  # nopep8

    _cli()
