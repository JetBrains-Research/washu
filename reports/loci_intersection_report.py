import os
import sys
import argparse
import datetime

from pathlib import Path
import pandas as pd


def _cli():
    # data_root = Path("/Volumes/BigData/bio")
    data_root = Path("/mnt/stripe/bio")
    loci_root = data_root / "raw-data/aging/loci_of_interest"
    # signal_root = data_root / "experiments/signal"

    ########################################################################
    parser = argparse.ArgumentParser(
        description="Generates intersection reports for predefined loci sets"
    )
    parser.add_argument('-o', '--out', required=True, metavar="PATH",
                        help="Output dir")
    parser.add_argument('-p', '--threads', type=int, default=4,
                        help="Threads number for parallel processing")
    parser.add_argument('--all', action="store_true",
                        help="Include outliers")
    parser.add_argument('--outliers', metavar="PATH",
                        default="/mnt/stripe/bio/experiments/aging/Y20O20.outliers.csv",
                        help="Outliers *.csv path")
    args = parser.parse_args()

    threads = args.threads
    outliers_df_path = args.outliers
    exclude_outliers = not args.all
    results_dir = Path(args.out)  # data_root / "experiments/aging/loci_of_interest.tables"
    ########################################################################

    loi_dict = loi.collect_loci(loci_root)

    golden_peaks_root = data_root / "experiments/aging/peak_calling"
    zinbra_peaks_root = data_root / "experiments/configs/Y20O20{}".format(
        "" if exclude_outliers else "_full"
    )

    peaks_map = {
        "zinbra": loi._collect_zinbra_peaks(zinbra_peaks_root),
        "golden": loi._collect_golden_peaks(golden_peaks_root, exclude_outliers)
    }

    outliers_df = None
    if outliers_df_path:
        outliers_df = pd.read_csv(outliers_df_path, delimiter="\t", skiprows=1,
                                  index_col="donor")
    ########################################################################

    print("----- [Report]: Default ----")
    report_default(loi_dict, results_dir, threads)
    print("----- [Report]: Consensus ----")
    report_consensus(loi_dict, results_dir, threads)
    print("----- [Report]: Repeats ----")
    report("repeats", loi_dict, results_dir, threads,
           consensus_yo=False, default=False)
    print("----- [Report]: Donors Zinbra ----")
    report_donors("zinbra", peaks_map, loi_dict, results_dir, threads, outliers_df)
    print("----- [Report]: Donors Golden ----")
    report_donors("golden", peaks_map, loi_dict, results_dir, threads, outliers_df)
    print("----- [Report]: Differential ChipSeq ----")
    report("chipseq_diff_loci", loi_dict, results_dir, threads)
    print("----- [Report]: Differential RnaSeq ----")
    report("rna_diff", loi_dict, results_dir, threads)
    print("----- [Report]: Pathways interesting ----")
    report("interesting_pathways", loi_dict, results_dir, threads,
           itself=False, chromhmm=False, default=False, repeats=False, consensus=False)

    # TODO: chromhmm -> fix files names
    # TODO: tfs


def report_default(loi_dict, outdir, threads):
    result_plot_path = outdir / "default.pdf"
    with PdfPages(str(result_plot_path)) as pdf:
        init_pdf_info(pdf)
        bm.process_intersection_metric(
            loi_dict['default'], loi_dict['default'],
            outdir / "default.csv", pdf,
            row_cluster=True, col_cluster=True, threads=threads, figsize=(20, 20),
            annotate_age=False, hist_mode=None)

        bm.process_intersection_metric(
            loi_dict['default'], loi_dict['chromhmm'],
            outdir / "default@chromhmm.csv", pdf,
            row_cluster=False, col_cluster=False, threads=threads, figsize=(12, 20),
            annotate_age=False, hist_mode=None)


def report_consensus(loi_dict, outdir, threads, consensus_type="median_consensus"):
    result_plot_path = outdir / "{}.pdf".format(consensus_type)
    with PdfPages(str(result_plot_path)) as pdf:
        init_pdf_info(pdf)
        consensus = loi_dict['zinbra_{}'.format(consensus_type)] \
            + loi_dict['golden_{}'.format(consensus_type)]
        bm.process_intersection_metric(
            consensus, loi_dict['default'],
            outdir / "{}@default.csv".format(consensus_type), pdf,
            row_cluster=True, col_cluster=True, threads=threads, figsize=(20, 15),
            annotate_age=False, hist_mode=None)

        bm.process_intersection_metric(
            consensus, consensus,
            outdir / "{}.csv".format(consensus_type), pdf,
            row_cluster=False, col_cluster=False, threads=threads, figsize=(12, 12),
            annotate_age=False, hist_mode=None)

        # YDS or ODS, but not "_ODS_without_YDS_median_consensus.bed"
        yo_consensus = [p for p in consensus if "DS" in p.name and "without" not in p.name]
        bm.process_intersection_metric(
            yo_consensus, yo_consensus,
            outdir / "{}_yo.csv".format(consensus_type), pdf,
            row_cluster=False, col_cluster=False, threads=threads, figsize=(10, 10),
            annotate_age=False, hist_mode=None)


def report_donors(tool, peaks_map, loi_dict, outdir, threads, outliers_df):
    peaks_dict = peaks_map[tool]

    result_plot_path = outdir / "{}_by_donor.pdf".format(tool)
    with PdfPages(str(result_plot_path)) as pdf:
        init_pdf_info(pdf)
        for hist in sorted(peaks_dict.keys()):
            bm.process_intersection_metric(
                peaks_dict[hist], loi_dict['default'],
                outdir / "{}_{}@default.csv".format(tool, hist), pdf,
                row_cluster=False, col_cluster=False, threads=threads, figsize=(20, 10),
                annotate_age=False, outliers_df=outliers_df, hist_mode=hist)


def report(key, loi_dict, outdir, threads, key_side_size=15, consensus_type="median_consensus",
           itself=True, chromhmm=True, default=True, repeats=True, consensus=True,
           consensus_yo=True):
    result_plot_path = outdir / "{}.pdf".format(key)

    with PdfPages(str(result_plot_path)) as pdf:
        init_pdf_info(pdf)
        if itself:
            bm.process_intersection_metric(
                loi_dict[key], loi_dict[key],
                outdir / "{}.csv".format(key), pdf,
                row_cluster=True, col_cluster=True, threads=threads,
                figsize=(key_side_size, key_side_size),
                annotate_age=False, hist_mode=None)

        # loci_key set:
        loci = []
        if chromhmm:
            loci.append(("chromhmm", 10))
        if default:
            loci.append(("default", 15))
        if repeats:
            loci.append(("repeats", 15))

        for loci_key, loci_side_size in loci:
            bm.process_intersection_metric(
                loi_dict[key], loi_dict[loci_key],
                outdir / "{}@{}.csv".format(key, loci_key), pdf,
                row_cluster=True, col_cluster=False, threads=threads,
                figsize=(loci_side_size, key_side_size),
                annotate_age=False, hist_mode=None)

            bm.process_intersection_metric(
                loi_dict[loci_key], loi_dict[key],
                outdir / "{}@{}.csv".format(loci_key, key), pdf,
                row_cluster=False, col_cluster=True, threads=threads,
                figsize=(key_side_size, loci_side_size),
                annotate_age=False, hist_mode=None)

        # consensus
        consensus_paths = loi_dict['zinbra_{}'.format(consensus_type)] \
            + loi_dict['golden_{}'.format(consensus_type)]
        if consensus:
            bm.process_intersection_metric(
                loi_dict[key], consensus_paths,
                outdir / "{}@{}.csv".format(key, consensus_type), pdf,
                row_cluster=True, col_cluster=False, threads=threads, figsize=(8, 15),
                annotate_age=False, hist_mode=None)

        # YDS or ODS, but not "_ODS_without_YDS_median_consensus.bed"
        if consensus_yo:
            consensus_yo_paths\
                = [p for p in consensus_paths if "DS" in p.name and "without" not in p.name]
            bm.process_intersection_metric(
                loi_dict[key], consensus_yo_paths,
                outdir / "{}@{}_yo.csv".format(key, consensus_type), pdf,
                row_cluster=True, col_cluster=False, threads=threads, figsize=(8, 15),
                annotate_age=False, hist_mode=None)


def init_pdf_info(pdf):
    # TODO: titles

    d = pdf.infodict()
    d['Title'] = 'Report: Intersection metric at different loci'
    d['Author'] = 'JetBrains Research BioLabs'
    d['Subject'] = 'outliers'
    # d['Keywords'] = 'outliers jetbrains aging'
    d['CreationDate'] = datetime.datetime.today()
    d['ModDate'] = datetime.datetime.today()


if __name__ == "__main__":
    parent_dir = os.path.dirname(os.path.realpath(__file__))
    project_root = os.path.abspath(os.path.join(parent_dir) + "/..")
    sys.path.insert(0, project_root)

    # Force matplotlib to not use any Xwindows backend.
    import matplotlib
    matplotlib.use('Agg')

    import reports.bed_metrics as bm
    import reports.loci_of_interest as loi

    from matplotlib.backends.backend_pdf import PdfPages

    _cli()
