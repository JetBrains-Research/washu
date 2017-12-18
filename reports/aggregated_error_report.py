import os
import sys
import re
import pandas as pd
from pathlib import Path

import argparse
import numpy as np
from itertools import chain


__author__ = 'petr.tsurinov@jetbrains.com'
help_data = """
Script creates pdf report with total error and peaks count plots from selected tools tuning 
parameters files for all histone modifications.
 1) Total error plot
 2) Peaks count bars
"""


def _cli():
    parser = argparse.ArgumentParser(description=help_data,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("folder", help="All histones folder")
    parser.add_argument("output", help="Output pdf path")
    parser.add_argument("tools", nargs='*', help="Tools folders")

    args = parser.parse_args()
    folder = Path(args.folder)
    pdf_path = args.output
    tools = args.tools

    with PdfPages(pdf_path) as pdf:
        plt.figure(figsize=(10, 15))
        threshold_map = {"H3K27ac": 0.3, "H3K27me3": 0.2, "H3K36me3": 0.45, "H3K4me3": 0.4,
                         "H3K4me1": 0.4}
        for hist_index, hist_mod in enumerate(["H3K27ac", "H3K27me3", "H3K36me3", "H3K4me3",
                                               "H3K4me1"]):
            plots = []
            current_tools = []
            colors = ['black', 'red', 'green', 'orange']
            ax = plt.subplot(5, 1, hist_index + 1)
            ax2 = ax.twinx()
            ax.set_ylabel(hist_mod + ' Total error')
            ax2.set_ylabel("Peaks count")
            ax2.grid(None)
            ax2.axes.get_xaxis().set_visible(False)

            for tool_index, tool in enumerate(tools):
                donor_peaks_map = peaks_count_map(folder / hist_mod / tool)
                parameters_path = folder / hist_mod / tool / "parameters.csv"

                if parameters_path.exists():
                    current_tools.append(tool)
                    parameters = pd.DataFrame.from_csv(path=str(parameters_path), index_col=None,
                                                       header=1, sep="\t")
                    parameters = parameters.drop_duplicates(['name', 'error'], keep='first')
                    n = len(parameters.index)
                    ind = np.arange(n)

                    plt.xticks(ind, parameters["name"].values)
                    plots.append(ax.plot(ind, parameters["error"], 'o-',
                                         color=colors[tool_index])[0])
                    ax2.bar(ind, [donor_peaks_map[donor] for donor in
                                  parameters["name"].values], 0.35, alpha=0.3,
                            color=colors[tool_index])

            plots.append(ax.plot(ind, [threshold_map[hist_mod]] * len(ind), ':', color='blue')[0])

            for label in ax.get_xmajorticklabels():
                label.set_rotation(90)

            plt.legend(plots, current_tools + ["Threshold"])

        plt.tight_layout()
        save_plot(pdf)


def peaks_count_map(path):
    peaks_paths = sorted(chain(path.glob("*Peak"), path.glob("*-island.bed"),
                               path.glob("*peaks.bed")))
    peaks_count = {}

    for peak_path in peaks_paths:
        match = re.match(".*([yo]d\d+).*", str(peak_path), flags=re.IGNORECASE)
        if match:
            peaks_count[match.group(1)] = Bed(str(peak_path)).count()
    return peaks_count


if __name__ == "__main__":
    # Force matplotlib to not use any Xwindows backend.
    import matplotlib

    matplotlib.use('Agg')

    parent_dir = os.path.dirname(os.path.realpath(__file__))
    project_root = os.path.abspath(os.path.join(parent_dir) + "/..")
    sys.path.insert(0, project_root)

    import matplotlib.pyplot as plt  # nopep8
    from matplotlib.backends.backend_pdf import PdfPages
    from reports.bed_metrics import save_plot  # nopep8
    from bed.bedtrace import Bed

    _cli()
