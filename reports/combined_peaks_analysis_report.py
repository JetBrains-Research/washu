import os
import sys
import re
import datetime
import pandas as pd
from pathlib import Path
from itertools import chain

# Force matplotlib to not use any Xwindows backend.
import matplotlib
matplotlib.use('Agg')

parent_dir = os.path.dirname(os.path.realpath(__file__))
project_root = os.path.abspath(os.path.join(parent_dir) + "/..")
sys.path.insert(0, project_root)

import seaborn as sns  # nopep8
from matplotlib.backends.backend_pdf import PdfPages  # nopep8
from bed.bedtrace import Bed  # nopep8
from reports.bed_metrics import color_annotator_chain, color_annotator_outlier, \
    color_annotator_age, bed_metric_table, plot_metric_heatmap, \
    label_converter_donor_and_tool  # nopep8
from reports.peak_metrics import calc_consensus, bar_consensus  # nopep8

__author__ = 'petr.tsurinov@jetbrains.com'
help_data = """
Usage: 
    combined_peaks_analysis_report.py [first peaks folder] [second peaks folder] [output pdf path]

Script creates pdf report with ChIP-seq peaks statistics:
 1) median peak consensus bar plot
 2) Metric #1 heatmap
"""
outliers_path = "/mnt/stripe/bio/experiments/aging/Y20O20.outliers.csv"
outliers_df = pd.read_csv(outliers_path, delimiter="\t", skiprows=1, index_col="donor")


def main():
    args = sys.argv

    if len(args) < 2:
        print(help_data)
        sys.exit(1)

    zinbra_folder_path = Path(args[1])
    standard_folder_path = Path(args[2])
    bed_files_paths = sorted(
        {str(zinbra_folder_path) + '/' + f for f in os.listdir(str(zinbra_folder_path)) if
         re.match('.*\.(?:broadPeak|bed|narrowPeak)$', f)})
    bed_files_paths.extend(
        {str(standard_folder_path) + '/' + f for f in os.listdir(str(standard_folder_path)) if
         re.match('.*\.(?:broadPeak|bed|narrowPeak)$', f)})
    bed_files_paths.sort()

    # tracks_paths = sorted(
    #     {bed_file for bed_file in bed_files_paths if re.match(".*([YO]D\d+).*", bed_file)})
    # od_paths_map = {re.findall('OD\\d+', track_path)[0] +
    #                 ("_zinbra" if str(zinbra_folder_path) in track_path else "_macs"):
    #                     Bed(track_path) for track_path in tracks_paths
    #                 if re.match('.*OD\\d+.*\.(?:broadPeak|bed|narrowPeak)$', track_path)}
    # yd_paths_map = {re.findall('YD\\d+', track_path)[0] +
    #                 ("_zinbra" if str(zinbra_folder_path) in track_path else "_macs"):
    #                     Bed(track_path) for track_path in tracks_paths
    #                 if re.match('.*YD\\d+.*\.(?:broadPeak|bed|narrowPeak)$', track_path)}

    anns = [color_annotator_age]
    hist_mod = re.match(".*(h3k\d{1,2}(?:me\d|ac)).*", str(zinbra_folder_path),
                        flags=re.IGNORECASE).group(1)
    if hist_mod in outliers_df.columns:
        anns.append(color_annotator_outlier(outliers_df, hist_mod))
    annotator = color_annotator_chain(*anns)

    peaks_paths = sorted(chain(zinbra_folder_path.glob("*golden*consensus*"),
                               zinbra_folder_path.glob("*zinbra*consensus*"),
                               zinbra_folder_path.glob("*Peak"),
                               zinbra_folder_path.glob("*-island.bed"),
                               zinbra_folder_path.glob("*peaks.bed")), key=loi.donor_order_id) + \
                  sorted(chain(standard_folder_path.glob("*golden*consensus*"),
                               standard_folder_path.glob("*zinbra*consensus*"),
                               standard_folder_path.glob("*Peak"),
                               standard_folder_path.glob("*-island.bed"),
                               standard_folder_path.glob("*peaks.bed")), key=loi.donor_order_id)
    df = bed_metric_table(peaks_paths, peaks_paths, threads=threads_num)
    df_no_out = df.copy()
    for donor in outliers_df.loc[:, hist_mod].index:
        if outliers_df.loc[:, hist_mod][donor] == 1:
            for df_index in df_no_out.index:
                if (donor + "_") in df_index or (donor + ".") in df_index:
                    del df_no_out[df_index]
                    df_no_out = df_no_out.drop(df_index)

    with PdfPages(args[3]) as pdf:
        # print("Calculating median consensus")
        # od_consensus_bed, yd_consensus_bed, yd_od_int_bed = \
        #     calc_consensus(od_paths_map, yd_paths_map, 2.0)
        # bar_consensus(od_paths_map, yd_paths_map, od_consensus_bed, yd_consensus_bed,
        #               yd_od_int_bed, threads_num, pdf)
        print("Calculating metric #1 indexes")
        sns.set(font_scale=0.5)
        for curr_df in (df, df_no_out):
            plot_metric_heatmap("Intersection metric", curr_df, figsize=(14, 14), save_to=pdf,
                                row_color_annotator=annotator, col_color_annotator=annotator,
                                row_label_converter=label_converter_donor_and_tool,
                                col_label_converter=label_converter_donor_and_tool)

        desc = pdf.infodict()
        desc['Title'] = 'Report: Combined peaks plots for different callers'
        desc['Author'] = 'JetBrains Research BioLabs'
        desc['Subject'] = 'peaks'
        desc['CreationDate'] = datetime.datetime.today()
        desc['ModDate'] = datetime.datetime.today()


if __name__ == "__main__":
    threads_num = 30

    import reports.loci_of_interest as loi

    main()
