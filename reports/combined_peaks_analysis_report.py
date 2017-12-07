import os
import sys
import re
import datetime
import argparse
import pandas as pd
from pathlib import Path

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


def _cli():
    parser = argparse.ArgumentParser(description=help_data,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("peaks_folder1", help="First tool peaks folder")
    parser.add_argument("peaks_folder2", help="Second tool peaks folder")
    parser.add_argument("output", help="Output pdf path")
    parser.add_argument('-p', '--threads', help="Threads number for parallel processing",
                        type=int, default=30)

    args = parser.parse_args()
    peaks_folder1 = Path(args.peaks_folder1)
    peaks_folder2 = Path(args.peaks_folder2)
    threads_num = args.threads
    pdf_path = args.output

    tracks_paths1 = sorted({peaks_folder1 / file for file in os.listdir(str(peaks_folder1)) if
                           re.match('.*(?:_peaks.bed|-island.bed|Peak)$', file)},
                           key=loi.donor_order_id)
    tracks_paths2 = sorted({peaks_folder2 / file for file in os.listdir(str(peaks_folder2)) if
                           re.match('.*(?:_peaks.bed|-island.bed|Peak)$', file)},
                           key=loi.donor_order_id)
    tracks_paths = tracks_paths1 + tracks_paths2
    tracks_names = list({str(tracks_path) for tracks_path in tracks_paths})

    od_paths_map = {re.findall('OD\\d+', track_name)[0] + _detect_tool(track_name):
                    track_name for track_name in tracks_names if re.match('.*OD\\d+.*', track_name)}
    yd_paths_map = {re.findall('YD\\d+', track_name)[0] + _detect_tool(track_name):
                    track_name for track_name in tracks_names if re.match('.*YD\\d+.*', track_name)}

    anns = [color_annotator_age]
    hist_mod = re.match(".*(h3k\d{1,2}(?:me\d|ac)).*", str(peaks_folder1),
                        flags=re.IGNORECASE).group(1)
    if hist_mod in outliers_df.columns:
        anns.append(color_annotator_outlier(outliers_df, hist_mod))
    annotator = color_annotator_chain(*anns)

    df = bed_metric_table(tracks_paths, tracks_paths, threads=threads_num)
    for donor in outliers_df.loc[:, hist_mod].index:
        if outliers_df.loc[:, hist_mod][donor] == 1:
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
        plot_metric_heatmap("Intersection metric", df, figsize=(14, 14), save_to=pdf,
                            row_color_annotator=annotator, col_color_annotator=annotator,
                            row_label_converter=label_converter_donor_and_tool,
                            col_label_converter=label_converter_donor_and_tool)

        desc = pdf.infodict()
        desc['Title'] = 'Report: Combined peaks plots for different callers'
        desc['Author'] = 'JetBrains Research BioLabs'
        desc['Subject'] = 'peaks'
        desc['CreationDate'] = datetime.datetime.today()
        desc['ModDate'] = datetime.datetime.today()


def _detect_tool(path):
    if "Peak" in path:
        return "_macs2"
    if "-island.bed" in path:
        return "_sicer"
    if "_peaks.bed" in path:
        return "_zinbra"
    return "_unknown"


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

    import reports.loci_of_interest as loi
    from matplotlib.backends.backend_pdf import PdfPages  # nopep8
    from reports.bed_metrics import color_annotator_chain, color_annotator_outlier, \
        color_annotator_age, bed_metric_table, plot_metric_heatmap, \
        label_converter_donor_and_tool  # nopep8
    from reports.peak_metrics import calc_consensus_file, bar_consensus  # nopep8

    _cli()
