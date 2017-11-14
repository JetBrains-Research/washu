import os
import sys
import re
import tempfile
import datetime
import pandas as pd
from pathlib import Path
from itertools import chain

# Force matplotlib to not use any Xwindows backend.
import matplotlib
matplotlib.use('Agg')

parent_dir = os.path.dirname(os.path.realpath(__file__))
project_root = os.path.abspath(os.path.join(parent_dir) + "/..")
sys.path.insert(0, project_root)

import seaborn as sns  # nopep8
from matplotlib.backends.backend_pdf import PdfPages  # nopep8
from scripts.util import regions_extension, age, is_od, is_yd, is_od_or_yd  # nopep8
from bed.bedtrace import Bed, run  # nopep8
from reports.bed_metrics import color_annotator_chain, color_annotator_outlier, \
    color_annotator_age, bed_metric_table, plot_metric_heatmap, \
    label_converter_donor_and_tool  # nopep8
from reports.peak_metrics import calc_consensus, venn_consensus, bar_consensus, \
    cumulative_consensus, calc_frip, frip_peaks, frip_boxplot, length_bar_plots  # nopep8

__author__ = 'petr.tsurinov@jetbrains.com'
help_data = """
Usage: peak_metrics.py [peaks folder] [output pdf path] [top peaks count (optional)]

Script creates pdf report with ChIP-seq peaks statistics:
 1) median peak consensus venn diagram
 2) median peak consensus bar plot
 3) cumulative consensus plot
 4) Jaccard similarity plot
 5) Jaccard index heatmap
 6) Frip/peaks plot
 7) Frip/age boxplot
 8) Peaks count/peaks length bar plots for each donor
"""
outliers_path = "/mnt/stripe/bio/experiments/aging/Y20O20.outliers.csv"
outliers_df = pd.read_csv(outliers_path, delimiter="\t", skiprows=1, index_col="donor")


def main():
    args = sys.argv

    if len(args) < 2:
        print(help_data)
        sys.exit(1)

    folder_path = Path(args[1])
    paths = sorted([str(f) for f in folder_path.iterdir() if regions_extension(f.name)])
    tmp_dir = Path(tempfile.gettempdir())
    filtered_paths = []

    if len(args) == 4:
        for path in paths:
            tmp_path = tmp_dir / "{}_{}.bed".format(Path(path).stem, args[3])
            with open(str(tmp_path), 'w') as f:
                run((["sort", "-k9nr", str(path)], ["head", "-n", args[3]]), stdout=f)
                filtered_paths.append(tmp_path.name)
    else:
        filtered_paths = paths

    tracks_paths = sorted({path for path in filtered_paths if is_od_or_yd(path)})
    od_paths_map = {age(track_path): Bed(track_path) for track_path in tracks_paths
                    if regions_extension(track_path) and is_od(track_path)}
    yd_paths_map = {age(track_path): Bed(track_path) for track_path in tracks_paths
                    if regions_extension(track_path) and is_yd(track_path)}
    rip_files = sorted([str(f) for f in folder_path.glob("*_rip.csv")])

    anns = [color_annotator_age]
    hist_mod = re.match(".*(h3k\d{1,2}(?:me\d|ac)).*", str(folder_path),
                        flags=re.IGNORECASE).group(1)
    if hist_mod in outliers_df.columns:
        anns.append(color_annotator_outlier(outliers_df, hist_mod))
    annotator = color_annotator_chain(*anns)

    peaks_paths = list(chain(folder_path.glob("*golden*consensus*"),
                             folder_path.glob("*zinbra*consensus*"), folder_path.glob("*Peak"),
                             folder_path.glob("*-island.bed"), folder_path.glob("*peaks.bed")))
    df = bed_metric_table(peaks_paths, peaks_paths, threads=threads_num)

    with PdfPages(args[2]) as pdf:
        print("Calculating median consensus")
        od_consensus_bed, yd_consensus_bed, yd_od_int_bed = \
            calc_consensus(od_paths_map, yd_paths_map, 2.0)
        venn_consensus(od_consensus_bed, yd_consensus_bed, 2.0, pdf)
        bar_consensus(od_paths_map, yd_paths_map, od_consensus_bed, yd_consensus_bed,
                      yd_od_int_bed, threads_num, pdf)
        print("Calculating cumulative consensus")
        cumulative_consensus(tracks_paths, pdf)
        print("Calculating jaccard indexes")
        sns.set(font_scale=0.75)
        plot_metric_heatmap("Intersection metric", df, figsize=(8, 8), save_to=pdf,
                            row_cluster=True, col_cluster=True,
                            row_color_annotator=annotator, col_color_annotator=annotator,
                            row_label_converter=label_converter_donor_and_tool,
                            col_label_converter=label_converter_donor_and_tool)
        print("Calculating frip vs age")
        age_labels, df = calc_frip(rip_files)
        frip_peaks(age_labels, df, pdf)
        frip_boxplot(age_labels, df, pdf)
        print("Calculating peaks count vs length")
        length_bar_plots(tracks_paths, 2.0, 4.0, threads_num, pdf)

        desc = pdf.infodict()
        desc['Title'] = 'Report: Peaks plots for data investigation'
        desc['Author'] = 'JetBrains Research BioLabs'
        desc['Subject'] = 'peaks'
        desc['CreationDate'] = datetime.datetime.today()
        desc['ModDate'] = datetime.datetime.today()


if __name__ == "__main__":
    threads_num = 30
    main()
