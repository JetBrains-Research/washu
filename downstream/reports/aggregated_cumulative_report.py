import argparse

import numpy as np
import pandas as pd

__author__ = 'petr.tsurinov@jetbrains.com'
help_data = """
Script creates pdf report with cumulative consensus plot for selected tools
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

    colors = ['black', 'red', 'green', 'orange']
    dfd = pd.read_csv(peaks_summary, sep='\t', comment='#')
    dfd = dfd[['donor', 'modification', 'tool', 'procedure', 'file', 'status']]

    with PdfPages(pdf_path) as pdf:
        for hist_index, hist_mod in enumerate(["H3K27ac", "H3K27me3", "H3K36me3", "H3K4me3",
                                               "H3K4me1"]):
            for tool_index, tool in enumerate(tools):
                plots = []
                curr_dfd = dfd.loc[np.logical_and(dfd['modification'] == hist_mod,
                                                  dfd['tool'] == tool)]
                if len(curr_dfd) > 0:
                    plt.figure()
                    for index, procedure in enumerate(["tuned", "default"]):
                        tracks_paths = curr_dfd.loc[
                            np.logical_and(curr_dfd['procedure'] == procedure,
                                           curr_dfd['status'] == 'ok')]['file']
                        if len(tracks_paths) > 0:
                            tracks_union = union(*[Bed(str(p)) for p in tracks_paths])
                            tracks_union.compute()

                            counts = [0] * len(tracks_paths)
                            for line in tracks_union.cat().split('\n'):
                                if line != '':
                                    parts = line.split("\t")
                                    count = len(parts[3].split("|"))
                                    counts[count - 1] += 1
                            counts.reverse()
                            counts_cumulative = list(np.cumsum(counts))
                            counts_cumulative.reverse()

                            plots.append(plt.plot(range(1, len(tracks_paths) + 1),
                                                  counts_cumulative, label="",
                                                  color=colors[index])[0])

                    plt.legend(plots, [tool + " tuned", tool + " default"])
                    plt.xlabel('Number of donors')
                    plt.ylabel('Number of peaks')
                    plt.title(hist_mod + " " + tool +
                              "\nreverse cumulative consensus peaks sum via number of donors")
                    plt.tight_layout()
                    pdf.savefig()
                    plt.close()


if __name__ == "__main__":
    # Force matplotlib to not use any Xwindows backend.
    import matplotlib
    matplotlib.use('Agg')

    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    from bed.bedtrace import Bed, union

    _cli()
