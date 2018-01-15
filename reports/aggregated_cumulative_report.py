import os
import sys
import argparse
import numpy as np
import pandas as pd
from pathlib import Path

__author__ = 'petr.tsurinov@jetbrains.com'
help_data = """
Script creates pdf report with cumulative consensus plot for selected tools
(all histone modifications).
"""
failed_tracks_path = "/mnt/stripe/bio/experiments/aging/Y20O20.failed_tracks.csv"
failed_tracks_df = pd.read_csv(failed_tracks_path, delimiter="\t", skiprows=1, index_col="donor")


def _cli():
    parser = argparse.ArgumentParser(description=help_data,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("folder_new", help="New results folder")
    parser.add_argument("folder_old", help="Old results folder")
    parser.add_argument("output", help="Output pdf path")
    parser.add_argument("tools", nargs='*', help="Tools folders")

    args = parser.parse_args()
    folder_new = Path(args.folder_new)
    folder_old = Path(args.folder_old)
    pdf_path = args.output
    tools = args.tools
    colors = ['black', 'red', 'green', 'orange']

    with PdfPages(pdf_path) as pdf:
        for hist_index, hist_mod in enumerate(["H3K27ac", "H3K27me3", "H3K36me3", "H3K4me3",
                                               "H3K4me1"]):
            for tool_index, tool in enumerate(tools):
                plots = []
                tool_new_path = folder_new / hist_mod / tool
                tool_old_path = folder_old / hist_mod / tool
                filter_donors = set(failed_tracks_df[failed_tracks_df[hist_mod] == 0].index.tolist())

                plt.figure()
                for index, tool_path in enumerate([tool_new_path, tool_old_path]):
                    tracks_paths = collect_peaks_in_folder(tool_path)
                    filter_paths = [p for p in tracks_paths if age(p.name) in filter_donors]
                    if len(filter_paths) > 0:
                        tracks_union = union(*[Bed(str(p)) for p in filter_paths])
                        tracks_union.compute()

                        counts = [0] * len(filter_paths)
                        for line in tracks_union.cat().split('\n'):
                            if line != '':
                                parts = line.split("\t")
                                count = len(parts[3].split("|"))
                                counts[count - 1] += 1
                        counts.reverse()
                        counts_cumulative = list(np.cumsum(counts))
                        counts_cumulative.reverse()

                        plots.append(plt.plot(range(1, len(filter_paths) + 1), counts_cumulative,
                                              label="", color=colors[index])[0])

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

    parent_dir = os.path.dirname(os.path.realpath(__file__))
    project_root = os.path.abspath(os.path.join(parent_dir) + "/..")
    sys.path.insert(0, project_root)

    import matplotlib.pyplot as plt  # nopep8
    from matplotlib.backends.backend_pdf import PdfPages
    from scripts.util import age
    from bed.bedtrace import Bed, union  # nopep8
    from reports.loci_of_interest import collect_peaks_in_folder

    _cli()
