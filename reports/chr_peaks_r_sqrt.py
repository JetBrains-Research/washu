import os
import re
import sys
import csv
import argparse
import numpy as np
from pathlib import Path
from itertools import chain
from sklearn import linear_model
from sklearn.metrics import r2_score

__author__ = 'petr.tsurinov@jetbrains.com'
help_data = """
Usage:
    chr_peaks_r_sqrt.py [input folder] [output folder]
    
    Calculates R^2 for linear regression constructed for peaks count on chr vs chr length.
    Plots graph with linear regression for each found peaks track.
"""
chr_sizes = "/mnt/stripe/bio/genomes/hg19/hg19.chrom.sizes"


def _cli():
    parser = argparse.ArgumentParser(description=help_data,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("peaks", help="Peaks folder")
    parser.add_argument("output", help="Output folder")

    args = parser.parse_args()
    input_path = Path(args.peaks)
    output_path = Path(args.output)

    peaks_paths = sorted(chain(input_path.glob("*Peak"),
                               input_path.glob("*-island.bed"),
                               input_path.glob("*peaks.bed")), key=loi.donor_order_id)

    print("Chr length vs peaks count")
    chr_lengths = {}
    peaks_count_template = {}
    with open(chr_sizes) as tsv:
        for line in csv.reader(tsv, dialect="excel-tab"):
            if re.match("^chr(\d{1,2}|[XY])$", line[0], flags=re.IGNORECASE):
                chr_lengths[line[0]] = int(line[1])
                peaks_count_template[line[0]] = 0

    for peaks_path in peaks_paths:
        peaks_count = peaks_count_template.copy()
        with open(peaks_path) as tsv:
            for line in csv.reader(tsv, dialect="excel-tab"):
                if line[0] in peaks_count.keys():
                    peaks_count[line[0]] = peaks_count[line[0]] + 1

        regr = linear_model.LinearRegression()
        regr.fit(np.asarray(list(chr_lengths.values())).reshape(-1, 1),
                 np.array(list(peaks_count.values())))
        peaks_count_pred = regr.predict(np.asarray(list(chr_lengths.values())).reshape(-1, 1))

        # Explained variance score: 1 is perfect prediction
        r2 = r2_score(list(peaks_count.values()), peaks_count_pred)
        print(peaks_path.stem + 'Variance score: %.2f' % r2)

        with PdfPages(str(output_path / ("%.2f_" % r2 + peaks_path.stem + ".pdf"))) as pdf:
            fig, ax = plt.subplots(figsize=(20, 5))
            ax.scatter(list(chr_lengths.values()), list(peaks_count.values()))
            ax.plot(list(chr_lengths.values()), peaks_count_pred, color='blue', linewidth=3)
            ax.legend(['Variance score: %.2f' % r2])
            ax.set_xlabel('Chromosome length', fontsize=18)
            ax.set_ylabel('Peaks count', fontsize=16)
            for i, txt in enumerate(chr_lengths.keys()):
                ax.annotate(txt, (list(chr_lengths.values())[i], list(peaks_count.values())[i]))
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
    import reports.loci_of_interest as loi
    from matplotlib.backends.backend_pdf import PdfPages

    _cli()
