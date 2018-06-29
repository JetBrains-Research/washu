import argparse

import numpy as np
import pandas as pd

__author__ = 'petr.tsurinov@jetbrains.com'
help_data = """
Script creates pdf report with bar plots for selected tools
(all histone modifications).
"""


def _cli():
    parser = argparse.ArgumentParser(description=help_data,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("peaks_summary", help="Peaks summary")
    parser.add_argument("output", help="Output pdf path")
    parser.add_argument("tools", nargs='*', help="Tools folders")

    args = parser.parse_args()
    peaks_summary = args.peaks_summary
    pdf_path = args.output
    tools = args.tools

    dfd = pd.read_csv(peaks_summary, sep='\t', comment='#')
    dfd = dfd[['modification', 'tool', 'file']].loc[dfd['procedure'] == 'tuned']

    with PdfPages(pdf_path) as pdf:
        for hist_index, hist_mod in enumerate(["H3K27ac", "H3K27me3", "H3K36me3", "H3K4me3",
                                               "H3K4me1"]):
            for tool_index, tool in enumerate(tools):
                plt.figure()
                curr_dfd = dfd.loc[np.logical_and(dfd['modification'] == hist_mod,
                                                  dfd['tool'] == tool)]

                if len(curr_dfd['file']) > 0:
                    bar_plot(hist_mod, tool, curr_dfd['file'], pdf, (10, 7), 14)


if __name__ == "__main__":
    # Force matplotlib to not use any Xwindows backend.
    import matplotlib
    matplotlib.use('Agg')

    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    from downstream.peak_metrics import bar_plot

    _cli()
