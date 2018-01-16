import argparse
import datetime
import os
import re
import sys
import tempfile
from itertools import chain
from pathlib import Path

import pandas as pd

__author__ = 'petr.tsurinov@jetbrains.com'
help_data = """
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
failed_tracks_path = "/mnt/stripe/bio/experiments/aging/Y20O20.failed_tracks.csv"
failed_tracks_df = pd.read_csv(failed_tracks_path, delimiter="\t", skiprows=1, index_col="donor")
loci_root = Path("/mnt/stripe/bio/raw-data/aging/loci_of_interest")


def _cli():
    parser = argparse.ArgumentParser(description=help_data,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("peaks", help="Peaks folder")
    parser.add_argument("output", help="Output pdf path")
    parser.add_argument("--count", type=int, help="Top peaks count")
    parser.add_argument('-p', '--threads', help="Threads number for parallel processing",
                        type=int, default=30)

    args = parser.parse_args()
    folder_path = Path(args.peaks)
    threads_num = args.threads
    pdf_path = args.output
    top_peaks_count = args.count

    paths = sorted([str(f) for f in folder_path.iterdir() if regions_extension(f.name)])
    tmp_dir = Path(tempfile.gettempdir())
    filtered_paths = []

    if top_peaks_count:
        for path in paths:
            tmp_path = tmp_dir / "{}_{}.bed".format(Path(path).stem, top_peaks_count)
            with open(str(tmp_path), 'w') as f:
                run((["sort", "-k9nr", str(path)], ["head", "-n", top_peaks_count]), stdout=f)
                filtered_paths.append(tmp_path.name)
    else:
        filtered_paths = paths

    tracks_paths = sorted({path for path in filtered_paths if is_od_or_yd(path)})
    od_paths_map = {age(track_path): track_path for track_path in tracks_paths
                    if regions_extension(track_path) and is_od(track_path)}
    yd_paths_map = {age(track_path): track_path for track_path in tracks_paths
                    if regions_extension(track_path) and is_yd(track_path)}
    rip_files = sorted([str(f) for f in folder_path.glob("*_rip.csv")])

    anns = [bm.color_annotator_age]
    hist_mod = re.match(".*(h3k\d{1,2}(?:me\d|ac)).*", str(folder_path),
                        flags=re.IGNORECASE).group(1)
    if hist_mod in failed_tracks_df.columns:
        anns.append(bm.color_annotator_outlier(failed_tracks_df, hist_mod))
    annotator = bm.color_annotator_chain(*anns)

    peaks_paths = sorted(chain(folder_path.glob("*sicer*consensus*"),
                               folder_path.glob("*macs2*consensus*"),
                               folder_path.glob("*zinbra*consensus*"), folder_path.glob("*Peak"),
                               folder_path.glob("*-island.bed"), folder_path.glob("*peaks.bed")),
                         key=loi.donor_order_id)

    df = bm.bed_metric_table(peaks_paths, peaks_paths, threads=threads_num)

    loci_dict = loi.collect_loci(loci_root)
    default_paths = []
    for key in ["top_level_paths", "enhancers", "regulatory", "repeats", 'chromhmm']:
        if key in loci_dict:
            default_paths.extend(loci_dict[key])
        else:
            print("Annotations not found:", str(loci_root / key), file=sys.stderr)
    loci_dict["default"] = sorted(default_paths, key=lambda p: p.name)

    df_loci = bm.bed_metric_table(peaks_paths, loci_dict['default'], threads=threads_num)

    with PdfPages(pdf_path) as pdf:
        print("Calculating median consensus")
        od_consensus_bed, yd_consensus_bed, yd_od_int_bed = \
            pm.calc_consensus_file(list(od_paths_map.values()), list(yd_paths_map.values()),
                                   percent=50)
        pm.venn_consensus(od_consensus_bed, yd_consensus_bed, 50, pdf)
        pm.bar_consensus(od_paths_map, yd_paths_map, od_consensus_bed, yd_consensus_bed,
                         yd_od_int_bed, threads_num, pdf)
        print("Calculating cumulative consensus")
        pm.cumulative_consensus(tracks_paths, pdf)
        print("Calculating Intersection metric")
        sns.set(font_scale=0.75)
        bm.plot_metric_heatmap("Intersection metric", df, figsize=(8, 8), save_to=pdf,
                               row_cluster=True, col_cluster=True,
                               row_color_annotator=annotator, col_color_annotator=annotator,
                               row_label_converter=bm.label_converter_donor_and_tool,
                               col_label_converter=bm.label_converter_donor_and_tool)
        bm.plot_metric_heatmap(
            "IM peaks@loci", df_loci, figsize=(15, 8), save_to=pdf,
            adjustments=dict(left=0.15, top=0.95, right=0.9, bottom=0.3),
            row_cluster=False, col_cluster=False,
            row_color_annotator=annotator,
            col_label_converter=loi.label_converter_shorten_loci,
            row_label_converter=bm.label_converter_donor_and_tool,
        )

        print("Calculating frip vs age")
        age_labels, df = pm.calc_frip(rip_files)
        pm.frip_peaks(age_labels, df, pdf)
        pm.frip_boxplot(age_labels, df, pdf)
        print("Calculating peaks count vs length")
        pm.length_bar_plots(tracks_paths, 2.0, 4.0, threads_num, pdf)

        desc = pdf.infodict()
        desc['Title'] = 'Report: Peaks plots for data investigation'
        desc['Author'] = 'JetBrains Research BioLabs'
        desc['Subject'] = 'peaks'
        desc['CreationDate'] = datetime.datetime.today()
        desc['ModDate'] = datetime.datetime.today()


if __name__ == "__main__":
    # Force matplotlib to not use any Xwindows backend.
    import matplotlib

    matplotlib.use('Agg')

    parent_dir = os.path.dirname(os.path.realpath(__file__))
    project_root = os.path.abspath(os.path.join(parent_dir) + "/..")
    sys.path.insert(0, project_root)

    import seaborn as sns  # nopep8
    from matplotlib.backends.backend_pdf import PdfPages
    from scripts.util import regions_extension, age, is_od, is_yd, is_od_or_yd
    from bed.bedtrace import run
    import downstream.bed_metrics as bm
    import downstream.loci_of_interest as loi
    import downstream.peak_metrics as pm

    _cli()
