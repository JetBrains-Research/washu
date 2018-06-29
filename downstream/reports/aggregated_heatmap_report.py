import argparse
import pandas as pd
import numpy as np
from pathlib import Path

__author__ = 'petr.tsurinov@jetbrains.com'
help_data = """
Script creates pdf report with heatmap consensus plot for selected tools
(all histone modifications).
"""


def _cli():
    parser = argparse.ArgumentParser(description=help_data,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("first_peaks_summary", help="First peaks summary")
    parser.add_argument("second_peaks_summary", help="Second peaks summary")
    parser.add_argument("output", help="Output pdf path")
    parser.add_argument("tools", nargs='*', help="Tools folders")
    parser.add_argument('-p', '--threads', help="Threads number for parallel processing",
                        type=int, default=6)

    args = parser.parse_args()
    first_peaks_summary = args.first_peaks_summary
    second_peaks_summary = args.second_peaks_summary
    pdf_path = args.output
    tools = args.tools
    threads_num = args.threads

    dfd = pd.read_csv(first_peaks_summary, sep='\t', comment='#')
    dfe = pd.read_csv(second_peaks_summary, sep='\t', comment='#')
    dfd = dfd[['donor', 'modification', 'tool', 'procedure', 'file', 'status']]
    dfe = dfe[['donor', 'modification', 'tool', 'procedure', 'file', 'status']]

    with PdfPages(pdf_path) as pdf:
        for hist_index, hist_mod in enumerate(["H3K27ac", "H3K27me3", "H3K36me3", "H3K4me3",
                                               "H3K4me1"]):
            for tool_index, tool in enumerate(tools):
                for procedure in ["default", "tuned"]:
                    curr_dfd = dfd.loc[np.logical_and(np.logical_and(dfd['tool'] == tool,
                                                      dfd['modification'] == hist_mod),
                                                      dfd['procedure'] == procedure)]
                    curr_dfe = dfe.loc[np.logical_and(np.logical_and(dfe['tool'] == tool,
                                                      dfe['modification'] == hist_mod),
                                                      dfe['procedure'] == procedure)]

                    plt.figure()
                    tracks_paths = [Path(f) for f in curr_dfd['file']] + \
                                   [Path(f) for f in curr_dfe['file']]
                    if len(tracks_paths) > 0:
                        df = bm.bed_metric_table(tracks_paths, tracks_paths, threads=threads_num)
                        failed_tracks_df = curr_dfd[['donor', 'status']].set_index('donor') \
                            .replace('ok', 0).replace('failed', 1)

                        anns = [bm.color_annotator_age,
                                bm.color_annotator_outlier(failed_tracks_df, 'status')]
                        annotator = bm.color_annotator_chain(*anns)

                        sns.set(font_scale=0.75)
                        bm.plot_metric_heatmap(
                            hist_mod + " " + tool + " " + procedure + " Overlap metric", df,
                            figsize=(8, 8), save_to=pdf,
                            row_cluster=True, col_cluster=True,
                            row_color_annotator=annotator,
                            col_color_annotator=annotator,
                            row_label_converter=bm.label_converter_donor_and_tool,
                            col_label_converter=bm.label_converter_donor_and_tool)


if __name__ == "__main__":
    # Force matplotlib to not use any Xwindows backend.
    import matplotlib

    matplotlib.use('Agg')

    import seaborn as sns
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    import downstream.bed_metrics as bm

    _cli()
