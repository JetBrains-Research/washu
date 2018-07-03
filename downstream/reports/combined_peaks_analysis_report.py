import argparse
import datetime
import re
from pathlib import Path

import pandas as pd
import numpy as np

__author__ = 'petr.tsurinov@jetbrains.com'
help_data = """
Usage:
    combined_peaks_analysis_report.py [first peaks summary] [second peaks summary] [output folder] 
    [first tool] [second tool]

Script creates pdf reports with ChIP-seq peaks statistics:
 1) median peak consensus bar plot
 2) overlap metric heatmap
"""


def _cli():
    parser = argparse.ArgumentParser(description=help_data,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("first_peaks_summary", help="First peaks summary")
    parser.add_argument("second_peaks_summary", help="Second peaks summary")
    parser.add_argument("output_folder", help="Output folder for pdf")
    parser.add_argument("tool1", help="First tool")
    parser.add_argument("tool2", help="Second tool")
    parser.add_argument('-p', '--threads', help="Threads number for parallel processing",
                        type=int, default=30)

    args = parser.parse_args()
    first_peaks_summary = args.first_peaks_summary
    second_peaks_summary = args.second_peaks_summary
    output_folder = Path(args.output_folder)
    tool1 = args.tool1
    tool2 = args.tool2
    threads_num = args.threads

    dfd = pd.read_csv(first_peaks_summary, sep='\t', comment='#')
    dfe = pd.read_csv(second_peaks_summary, sep='\t', comment='#')
    dfd = dfd[['donor', 'modification', 'tool', 'procedure', 'file', 'status']]
    dfe = dfe[['donor', 'modification', 'tool', 'procedure', 'file', 'status']]

    today = datetime.datetime.today()
    pdf_path = str(output_folder / (tool1 + "_vs_" + tool2 + "_" +
                                    today.strftime("%d.%m.%Y") + ".pdf"))
    with PdfPages(pdf_path) as pdf:
        for hist_mod in ["H3K27ac", "H3K27me3", "H3K36me3", "H3K4me3", "H3K4me1"]:
            for procedure in ["default", "tuned"]:
                curr_dfd1 = dfd.loc[
                    np.logical_and(np.logical_and(np.logical_and(dfd['tool'] == tool1,
                                                                 dfd['modification'] == hist_mod),
                                                  dfd['procedure'] == procedure),
                                   dfd['status'] == "ok")]
                curr_dfe1 = dfe.loc[
                    np.logical_and(np.logical_and(np.logical_and(dfe['tool'] == tool1,
                                                                 dfe['modification'] == hist_mod),
                                                  dfe['procedure'] == procedure),
                                   dfe['status'] == "ok")]
                curr_dfd2 = dfd.loc[
                    np.logical_and(np.logical_and(np.logical_and(dfd['tool'] == tool2,
                                                                 dfd['modification'] == hist_mod),
                                                  dfd['procedure'] == procedure),
                                   dfd['status'] == "ok")]
                curr_dfe2 = dfe.loc[
                    np.logical_and(np.logical_and(np.logical_and(dfe['tool'] == tool2,
                                                                 dfe['modification'] == hist_mod),
                                                  dfe['procedure'] == procedure),
                                   dfe['status'] == "ok")]

                tracks_paths1 = sorted([Path(f) for f in curr_dfd1['file']] +
                                       [Path(f) for f in curr_dfe1['file']], key=loi.donor_order_id)
                tracks_paths2 = sorted([Path(f) for f in curr_dfd2['file']] +
                                       [Path(f) for f in curr_dfe2['file']], key=loi.donor_order_id)
                if len(tracks_paths1) > 0 and len(tracks_paths2) > 0:
                    tracks_paths = tracks_paths1 + tracks_paths2

                    tracks_names = list({str(tracks_path) for tracks_path in tracks_paths})
                    od_paths_map = {re.findall('OD\\d+', track_name)[0] + detect_tool(track_name):
                                    track_name for track_name in tracks_names if
                                    re.match('.*OD\\d+.*', track_name)}
                    yd_paths_map = {re.findall('YD\\d+', track_name)[0] + detect_tool(track_name):
                                    track_name for track_name in tracks_names if
                                    re.match('.*YD\\d+.*', track_name)}
                    df = bed_metric_table(tracks_paths, tracks_paths, threads=threads_num)

                    print("Calculating median consensus")
                    od_consensus_bed, yd_consensus_bed, yd_od_int_bed = \
                        calc_consensus_file(list(od_paths_map.values()),
                                            list(yd_paths_map.values()), percent=50)
                    bar_consensus(od_paths_map, yd_paths_map, od_consensus_bed, yd_consensus_bed,
                                  yd_od_int_bed, threads_num, pdf, (10, 10), 10)

                    print("Calculating metric #1 indexes")
                    sns.set(font_scale=0.75)
                    g = plot_metric_heatmap(
                        hist_mod + " " + tool1 + " vs " + tool2 + " " + procedure +
                        " Overlap metric", df, figsize=(14, 14),
                        save_to=pdf,
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


if __name__ == "__main__":
    # Force matplotlib to not use any Xwindows backend.
    import matplotlib

    matplotlib.use('Agg')

    import seaborn as sns
    import downstream.loci_of_interest as loi
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from downstream.bed_metrics import bed_metric_table, plot_metric_heatmap, \
        label_converter_donor_and_tool, save_plot
    from downstream.peak_metrics import calc_consensus_file, bar_consensus, detect_tool

    _cli()
