import argparse
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from downstream.bed_metrics import save_plot  # nopep8

__author__ = 'petr.tsurinov@jetbrains.com'
help_data = """
Script creates pdf report with error plots from two selected tools tuning parameters files:
 1) Total error
 2) Peak start error
 3) Peak end error
 4) No peaks error
 5) Peaks error
"""


def _cli():
    parser = argparse.ArgumentParser(description=help_data,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("parameters1", help="First tool parameters file")
    parser.add_argument("parameters2", help="Second tool parameters file")
    parser.add_argument("output", help="Output pdf path")

    args = parser.parse_args()
    parameters_path1 = Path(args.parameters1)
    parameters_path2 = Path(args.parameters2)
    pdf_path = args.output

    parameters1 = pd.DataFrame.from_csv(path=str(parameters_path1), index_col=None, header=1,
                                        sep="\t")
    parameters1 = parameters1.drop_duplicates(['name', 'error'], keep='first')
    parameters2 = pd.DataFrame.from_csv(path=str(parameters_path2), index_col=None, header=1,
                                        sep="\t")
    parameters2 = parameters2.drop_duplicates(['name', 'error'], keep='first')
    n = max(len(parameters1.index), len(parameters2.index))
    ind = np.arange(n)

    with PdfPages(pdf_path) as pdf:
        plt.figure(figsize=(10, 15))

        for index, (error, label) in enumerate([("error", 'Total error'),
                                                ("error_S", 'Peak start error'),
                                                ("error_E", 'Peak end error'),
                                                ("error_N", 'No peaks error'),
                                                ("error_P", 'Peaks error')]):
            plt.subplot(5, 1, index + 1)
            error_plot(ind, parameters_path1, parameters1[error], parameters_path2,
                       parameters2[error], parameters1["name"].values, label)

        plt.tight_layout()
        save_plot(pdf)


def error_plot(ind, parameters_path1, errors1, parameters_path2, errors2, names, label):
    p1 = plt.plot(ind, errors1, 'o-', color='black')
    p2 = plt.plot(ind, errors2, 'o-', color='red')
    plt.xticks(ind, names, rotation=90)
    plt.ylabel(label)
    plt.legend((p1[0], p2[0]), (detect_tool(parameters_path1), detect_tool(parameters_path2)))


def detect_tool(path):
    return re.search('(zinbra|macs_broad|macs_narrow|sicer)', str(path),
                     flags=re.IGNORECASE).group(0)


if __name__ == "__main__":
    # Force matplotlib to not use any Xwindows backend.
    import matplotlib

    matplotlib.use('Agg')

    parent_dir = os.path.dirname(os.path.realpath(__file__))
    project_root = os.path.abspath(os.path.join(parent_dir) + "/..")
    sys.path.insert(0, project_root)

    import matplotlib.pyplot as plt  # nopep8
    from matplotlib.backends.backend_pdf import PdfPages

    _cli()
