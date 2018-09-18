import argparse
import re
from itertools import chain
from pathlib import Path

import numpy as np
import pandas as pd

__author__ = 'petr.tsurinov@jetbrains.com'
help_data = """
Script creates pdf report with total error, peaks count and frip plots from selected tools tuning
parameters files for all histone modifications.
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
        fig1 = plt.figure(figsize=(10, 15))
        fig2 = plt.figure(figsize=(10, 15))
        threshold_map = {"H3K27ac": 0.28, "H3K27me3": 0.1, "H3K36me3": 0.24, "H3K4me3": 0.13,
                         "H3K4me1": 0.28}

        for hist_index, hist_mod in enumerate(["H3K27ac", "H3K27me3", "H3K36me3", "H3K4me3",
                                               "H3K4me1"]):
            max_n = 0
            current_tools1 = {}
            current_tools2 = {}
            colors = ['black', 'red', 'green', 'orange']

            ax1_1 = fig1.add_subplot(5, 1, hist_index + 1)
            ax1_1.set_ylabel(hist_mod + ' Total error')

            ax1_2 = ax1_1.twinx()
            ax1_2.set_ylabel("Peaks count")
            ax1_2.grid(None)
            ax1_2.axes.get_xaxis().set_visible(False)

            ax2_1 = fig2.add_subplot(5, 1, hist_index + 1)
            ax2_1.set_ylabel(hist_mod + ' FRiP')

            ax2_2 = ax2_1.twinx()
            ax2_2.set_ylabel("Peaks count")
            ax2_2.grid(None)
            ax2_2.axes.get_xaxis().set_visible(False)

            for tool_index, tool in enumerate(tools):
                curr_dfd = dfd.loc[np.logical_and(dfd['modification'] == hist_mod,
                                                  dfd['tool'] == tool)]
                if len(curr_dfd) < 1:
                    continue
                tool_path = Path(curr_dfd['file'].iloc[0].rsplit('/', 1)[0])
                tool_name = str(tool_path).rsplit('/', 1)[1]
                parameters_path = tool_path / "{}_{}_parameters.csv".format(hist_mod, tool_name)
                if parameters_path.exists():
                    donor_peaks_map = peaks_count_map(tool_path)

                    parameters = pd.read_csv(
                        path=str(parameters_path),
                        index_col=None, header=1, sep="\t"
                    )
                    parameters = parameters.drop_duplicates(['name', 'error'], keep='first')
                    n = len(parameters.index)
                    ind = np.arange(n)

                    ax1_1.set_xticks(ind)
                    ax1_1.set_xticklabels(parameters["name"].values)
                    current_tools1[tool] = ax1_1.plot(ind, parameters["error"], 'o-',
                                                      color=colors[tool_index])[0]
                    ax1_2.bar(ind, [donor_peaks_map[donor] for donor in
                                    parameters["name"].values], 0.35, alpha=0.3,
                              color=colors[tool_index])

                    rip_files = sorted([str(f) for f in tool_path.glob("*_rip.csv")])
                    if len(rip_files) > 0:
                        age, frip_df = pm.calc_frip(rip_files)

                        n = len(frip_df.index)
                        max_n = max(max_n, n)
                        ind = np.arange(n)

                        names = sorted(frip_df.index, key=lambda name: (name[:2], int(name[2:])))
                        ax2_1.set_xticks(ind)
                        ax2_1.set_xticklabels(names)
                        frip_df_zeros = frip_df["frip"][names]
                        frip_df_zeros[np.isnan(frip_df_zeros)] = 0
                        current_tools2[tool] = ax2_1.plot(ind, frip_df_zeros, 'o-',
                                                          color=colors[tool_index])[0]
                        ax2_2.bar(ind, [donor_peaks_map[donor] for donor in names], 0.35,
                                  alpha=0.3, color=colors[tool_index])

            ind = np.arange(max_n)
            current_tools1["Threshold"] = ax1_1.plot(ind, [threshold_map[hist_mod]] * len(ind),
                                                     ':', color='blue')[0]

            for label in ax1_1.get_xmajorticklabels():
                label.set_rotation(90)
            for label in ax2_1.get_xmajorticklabels():
                label.set_rotation(90)

            ax1_1.legend(current_tools1.values(), current_tools1.keys())
            ax2_1.legend(current_tools2.values(), current_tools2.keys())

        fig1.tight_layout()
        fig2.tight_layout()
        pdf.savefig(fig1)
        pdf.savefig(fig2)
        plt.close()


def peaks_count_map(path):
    peaks_paths = sorted(chain(path.glob("*Peak"), path.glob("*-island.bed"),
                               path.glob("*peaks.bed")))
    peaks_count = {}

    for peak_path in peaks_paths:
        match = re.match(r".*([yo]d\d+).*", str(peak_path), flags=re.IGNORECASE)
        if match:
            peaks_count[match.group(1)] = Bed(str(peak_path)).count()
    return peaks_count


if __name__ == "__main__":
    # Force matplotlib to not use any Xwindows backend.
    import matplotlib
    matplotlib.use('Agg')

    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    import downstream.peak_metrics as pm
    from bed.bedtrace import Bed

    _cli()
