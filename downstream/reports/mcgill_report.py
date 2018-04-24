import argparse
import datetime
import re
import sys
import tempfile
from itertools import chain
from pathlib import Path
import matplotlib.pyplot as plt
from downstream.bed_metrics import save_plot
import multiprocessing
import functools
from collections import namedtuple

import pandas as pd

__author__ = 'alexey.dievsky@jetbrains.com'
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

loci_root = Path("/mnt/stripe/bio/raw-data/aging/loci_of_interest")

Group = namedtuple('Group', 'name color')
MONO = Group('M', 'blue')
TCELL = Group('T', 'red')
UNKNOWN = Group('', 'black')

def mcgill_group(c):
    if re.match('.*mn\\d+.*', str(c), flags=re.IGNORECASE):
        return MONO
    if re.match('.*tc\\d+.*', str(c), flags=re.IGNORECASE):
        return TCELL
    return UNKNOWN

def is_mono_or_tcell(c):
    return is_mono(c) or is_tcell(c)

def is_mono(c):
    return mcgill_group(c) == MONO

def is_tcell(c):
    return mcgill_group(c) == TCELL

def calc_consensus_file(mono_files_paths, tcell_files_paths, count=0, percent=0):
    mono_cons = consensus(mono_files_paths, count, percent)
    tcell_cons = consensus(tcell_files_paths, count, percent)

    mono_consensus_path = tempfile.NamedTemporaryFile(mode='wb', suffix='.bed', prefix='mn_consensus',
                                                    delete=False)
    tcell_consensus_path = tempfile.NamedTemporaryFile(mode='wb', suffix='.bed', prefix='tc_consensus',
                                                    delete=False)
    mono_consensus_path.write(mono_cons)
    tcell_consensus_path.write(tcell_cons)
    mono_consensus_path.close()
    tcell_consensus_path.close()

    mono_consensus_bed = Bed(mono_consensus_path.name)
    tcell_consensus_bed = Bed(tcell_consensus_path.name)
    mono_tcell_int_bed = intersect(mono_consensus_bed, tcell_consensus_bed)
    mono_tcell_int_bed.compute()

    return mono_consensus_bed, tcell_consensus_bed, mono_tcell_int_bed

def venn_consensus(mono_consensus_bed, tcell_consensus_bed, percent, save_to=None):
    """
    Plots venn diagram for consensus of selected scale:

    :param od_consensus_bed: OD bed with od group consensus
    :param yd_consensus_bed: YD bed with yd group consensus
    :param percent: percent of tracks count needed for consensus
    :param save_to: Object for plots saving
    """
    plt.figure()
    plt.title("Required consensus: %.2f%%" % percent)
    metapeaks({'Monocytes': mono_consensus_bed, 'T cells': tcell_consensus_bed})
    save_plot(save_to)

def donor_order_id(donor_data):
    chunks = donor_data[0].split('_')
    donor_id = chunks[0]
    if len(chunks) > 1:
        return chunks[1], donor_id[:2], int(donor_id[2:])
    return donor_id

def bar_consensus(mono_paths_map, tcell_paths_map, mono_consensus_bed, tcell_consensus_bed, mono_tcell_int_bed,
                  threads_num, save_to=None, figsize=(10, 10), fontsize=6):
    """
    Plots venn diagram and bar plot for consensus of selected scale:

    :param mono_paths_map: OD names as keys, od beds as values
    :param tcell_paths_map: YD names as keys, yd beds as values
    :param mono_consensus_bed: OD bed with od group consensus
    :param tcell_consensus_bed: YD bed with yd group consensus
    :param mono_tcell_int_bed: BED with intersection of od and yd groups consensuses
    :param figsize: Plot figure size
    :param save_to: Object for plots saving
    :param fontsize: Size of xlabels on plot
    :param threads_num: Threads number for parallel execution
    """
    pool = multiprocessing.Pool(processes=threads_num)
    n = len(tcell_paths_map) + len(mono_paths_map)
    ind = np.arange(n)

    tcell_result = pool.map(
        functools.partial(pm.groups_sizes, common_bed=mono_tcell_int_bed, own_group_bed=tcell_consensus_bed,
                          opposite_group_bed=mono_consensus_bed), sorted(tcell_paths_map.items(),
                                                                         key=operator.itemgetter(0)))
    mono_result = pool.map(
        functools.partial(pm.groups_sizes, common_bed=mono_tcell_int_bed, own_group_bed=mono_consensus_bed,
                          opposite_group_bed=tcell_consensus_bed), sorted(mono_paths_map.items(),
                                                                          key=operator.itemgetter(0)))
    result = sorted(tcell_result + mono_result, key=donor_order_id)
    result_columns = list(zip(*result))

    plt.figure(figsize=figsize)
    width = 0.35
    p1 = plt.bar(ind, result_columns[1], width, color='green')
    p2 = plt.bar(ind, result_columns[2], width, bottom=[mono_tcell_int_bed.count()] * n, color='blue')
    p3 = plt.bar(ind, result_columns[3], width, bottom=[mono_tcell_int_bed.count() +
                                                        max(result_columns[2])] * n, color='orange')
    p4 = plt.bar(ind, result_columns[4], width, bottom=[mono_tcell_int_bed.count() +
                                                        max(result_columns[2]) +
                                                        max(result_columns[3])] * n, color='black')
    plt.ylabel('Peaks count')
    plt.xticks(ind, result_columns[0], rotation=90, fontsize=fontsize)
    plt.legend((p1[0], p2[0], p3[0], p4[0]),
               ('Common', 'Own Group', 'Opposite Group', 'Individual'))
    plt.tight_layout()
    save_plot(save_to)

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

    tracks_paths = sorted({path for path in filtered_paths if is_mono_or_tcell(path)})
    mono_paths_map = {donor(track_path): track_path for track_path in tracks_paths
                    if regions_extension(track_path) and is_mono(track_path)}
    tcell_paths_map = {donor(track_path): track_path for track_path in tracks_paths
                    if regions_extension(track_path) and is_tcell(track_path)}
    rip_files = sorted([str(f) for f in folder_path.glob("*_rip.csv")])

    anns = [bm.color_annotator_age]
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
        mono_consensus_bed, tcell_consensus_bed, mono_tcell_int_bed = \
            calc_consensus_file(list(mono_paths_map.values()), list(tcell_paths_map.values()),
                                   percent=50)
        venn_consensus(mono_consensus_bed, tcell_consensus_bed, 50, pdf)
        bar_consensus(mono_paths_map, tcell_paths_map, mono_consensus_bed, tcell_consensus_bed,
                         mono_tcell_int_bed, threads_num, pdf)
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

    import seaborn as sns
    from matplotlib.backends.backend_pdf import PdfPages
    from downstream.aging import regions_extension, donor
    from bed.bedtrace import run, consensus, metapeaks, intersect
    import downstream.bed_metrics as bm
    import downstream.loci_of_interest as loi
    import downstream.peak_metrics as pm

    _cli()