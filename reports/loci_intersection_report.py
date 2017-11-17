import os
import sys
import argparse
import datetime

from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


def _cli():
    # data_root = Path("/Volumes/BigData/bio")
    data_root = Path("/mnt/stripe/bio")
    loci_root = data_root / "raw-data/aging/loci_of_interest"
    # signal_root = data_root / "experiments/signal"

    ########################################################################
    parser = argparse.ArgumentParser(
        description="Generates intersection reports for predefined loci sets",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
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
    results_dir.mkdir(parents=True, exist_ok=True)
    ########################################################################

    loci_dict = loi.collect_loci(loci_root)

    # default peaks
    golden_peaks_root = data_root / "experiments/aging/peak_calling"
    zinbra_peaks_root = data_root / "experiments/configs/Y20O20{}/peaks".format(
        "" if exclude_outliers else "_full"
    )
    # tuned peaks
    zinbra_peaks_root_tuned = data_root / "experiments/configs/Y20O20_full/benchmark_peaks"
    golden_peaks_root_tuned = data_root / "experiments/aging/peak_calling/benchmark_peaks"

    peaks_map = {
        "zinbra": loi._collect_zinbra_peaks(zinbra_peaks_root),
        "golden": loi._collect_golden_peaks(golden_peaks_root, exclude_outliers),
        "zinbra_tuned": loi._collect_zinbra_peaks(zinbra_peaks_root_tuned),
        "golden_tuned": loi._collect_golden_peaks(golden_peaks_root_tuned, None)
    }

    outliers_df = None
    if outliers_df_path:
        outliers_df = pd.read_csv(outliers_df_path, delimiter="\t", skiprows=1,
                                  index_col="donor")
    ########################################################################

    print("----- [Report]: Default ----")
    report_default(loci_dict, results_dir, threads)
    print("----- [Report]: Consensus ----")
    report_consensus(loci_dict, results_dir, threads)

    for key in sorted(peaks_map.keys()):
        print("----- [Report]: Donors {} ----".format(key))
        report_donors(key, peaks_map, loci_dict, results_dir, threads, outliers_df)

    print("----- [Report]: Repeats ----")
    report("repeats", loci_dict, results_dir, threads,
           consensus_yo=False, default=False)
    print("----- [Report]: TFs ----")
    report("tfs", loci_dict, results_dir, threads,
           consensus_yo=False, default=False)
    print("----- [Report]: Differential ChipSeq ----")
    report("chipseq_diff_loci", loci_dict, results_dir, threads)
    print("----- [Report]: Differential RnaSeq ----")
    report("rna_diff", loci_dict, results_dir, threads)

    print("----- [Report]: Pathways interesting ----")
    report("interesting_pathways", loci_dict, results_dir, threads,
           key_side_size=200,
           itself=False, chromhmm=False, default=False, repeats=False, consensus=False)

    print("----- [Report]: Pathways NOTCH ----")
    notch_pathways = [p for p in (loci_root / "interesting_pathways").glob('R-HSA-266082*.bed')]
    notch_pathways.extend(
        [p for p in (loci_root / "interesting_pathways").glob('R-HSA-1912399*.bed')])
    notch_pathways.extend(
        [p for p in (loci_root / "interesting_pathways").glob('R-HSA-264460*_cds.bed')])
    notch_pathways.sort(key=lambda p: p.name)
    loci_dict["notch_pathways"] = notch_pathways

    report("notch_pathways", loci_dict, results_dir, threads,
           key_side_size=20,
           itself=False, chromhmm=False, default=False, repeats=False, consensus=False)

    # Stat tests:
    loci = set(loci_dict.keys())
    # lets check large pathways and all data after all other loci:
    loci.remove(None)
    loci.remove("other_pathways")
    loci.remove("interesting_pathways")
    loci = sorted(loci)
    # loci.extend(["interesting_pathways", "other_pathways", None])
    loci.extend(["interesting_pathways", "other_pathways"])
    for loci_key in loci:
        for key in ["zinbra_tuned", "golden_tuned"]:
            print("----- [Stat tests]: Donors {}@{} ----".format(key, loci_key))
            test_donors(key, peaks_map, loci_dict, loci_key, results_dir, threads, outliers_df)


def _adjustment():
    return dict(left=0.15, top=0.95, right=0.65, bottom=0.25)


def _adjustment_wrc():
    return dict(left=0.15, top=0.95, right=0.62, bottom=0.3)


def report_default(loci_dict, outdir, threads):
    result_plot_path = outdir / "plot_default.pdf"
    with PdfPages(str(result_plot_path)) as pdf:
        init_pdf_info(pdf)
        process_intersection_metric(
            loci_dict['default'], loci_dict['default'],
            outdir / "default.csv", pdf,
            adjustments=_adjustment(),
            col_label_converter=loi.label_converter_shorten_loci,
            row_label_converter=loi.label_converter_shorten_loci,
            row_cluster=False, col_cluster=False, threads=threads, figsize=(20, 16))

        process_intersection_metric(
            loci_dict['default'], loci_dict['chromhmm'],
            outdir / "default@chromhmm.csv", pdf,
            adjustments=_adjustment(),
            col_label_converter=loi.label_converter_shorten_loci,
            row_label_converter=loi.label_converter_shorten_loci,
            row_cluster=False, col_cluster=False, threads=threads, figsize=(15, 15))


def report_consensus(loci_dict, outdir, threads, consensus_type="median_consensus"):
    result_plot_path = outdir / "plot_{}.pdf".format(consensus_type)
    with PdfPages(str(result_plot_path)) as pdf:
        init_pdf_info(pdf)
        consensus = sorted(
            loci_dict['zinbra_{}'.format(consensus_type)] + loci_dict['golden_{}'.format(
                consensus_type)],
            key=lambda f: f.name
        )
        process_intersection_metric(
            consensus, loci_dict['default'],
            outdir / "{}@default.csv".format(consensus_type), pdf,
            col_label_converter=loi.label_converter_shorten_loci,
            row_label_converter=loi.label_converter_shorten_loci,
            adjustments=_adjustment_wrc(),
            row_cluster=False, col_cluster=True, threads=threads, figsize=(20, 15))

        process_intersection_metric(
            consensus, loci_dict['chromhmm'],
            outdir / "{}@chromhmm.csv".format(consensus_type), pdf,
            adjustments=_adjustment(),
            col_label_converter=loi.label_converter_shorten_loci,
            row_label_converter=loi.label_converter_shorten_loci,
            row_cluster=False, col_cluster=False, threads=threads, figsize=(20, 15))

        process_intersection_metric(
            consensus, consensus,
            outdir / "{}.csv".format(consensus_type), pdf,
            adjustments=_adjustment(),
            col_label_converter=loi.label_converter_shorten_loci,
            row_label_converter=loi.label_converter_shorten_loci,
            row_cluster=False, col_cluster=False, threads=threads, figsize=(20, 15))

        process_intersection_metric(
            consensus, loci_dict['repeats'],
            outdir / "{}@repeats.csv".format(consensus_type), pdf,
            col_label_converter=loi.label_converter_shorten_loci,
            row_label_converter=loi.label_converter_shorten_loci,
            adjustments=_adjustment_wrc(),
            row_cluster=False, col_cluster=True, threads=threads, figsize=(20, 15))

        # YDS or ODS, but not "_ODS_without_YDS_median_consensus.bed"
        yo_consensus = [p for p in consensus if "DS" in p.name and "without" not in p.name]
        process_intersection_metric(
            yo_consensus, yo_consensus,
            outdir / "{}_yo.csv".format(consensus_type), pdf,
            adjustments=dict(left=0.15, top=0.95, right=0.75, bottom=0.25),
            col_label_converter=loi.label_converter_shorten_loci,
            row_label_converter=loi.label_converter_shorten_loci,
            row_cluster=False, col_cluster=False, threads=threads, figsize=(10, 10))


def report_donors(tool, peaks_map, loci_dict, outdir, threads, outliers_df):
    peaks_dict = peaks_map[tool]

    result_plot_path = outdir / "plot_{}_by_donor.pdf".format(tool)
    with PdfPages(str(result_plot_path)) as pdf:
        init_pdf_info(pdf)
        for hist in sorted(peaks_dict.keys()):
            anns = [bm.color_annotator_age]
            if hist and outliers_df is not None:
                if hist in outliers_df.columns:
                    anns.append(bm.color_annotator_outlier(outliers_df, hist))
            annotator = None if not anns else bm.color_annotator_chain(*anns)

            process_intersection_metric(
                peaks_dict[hist], loci_dict['default'],
                outdir / "{}_{}@default.csv".format(tool, hist), pdf,
                adjustments=dict(left=0.15, top=0.95, right=0.9, bottom=0.3),
                row_cluster=False, col_cluster=False, threads=threads, figsize=(20, 15),
                row_color_annotator=annotator,
                col_label_converter=loi.label_converter_shorten_loci,
                row_label_converter=bm.label_converter_donor_and_tool,
            )


def report(key, loci_dict, outdir, threads, key_side_size=15, consensus_type="median_consensus",
           itself=True, chromhmm=True, default=True, repeats=True, consensus=True,
           consensus_yo=True):
    result_plot_path = outdir / "plot_{}.pdf".format(key)

    with PdfPages(str(result_plot_path)) as pdf:
        init_pdf_info(pdf)
        if itself:
            process_intersection_metric(
                loci_dict[key], loci_dict[key],
                outdir / "{}.csv".format(key), pdf,
                adjustments=_adjustment_wrc(),
                col_label_converter=loi.label_converter_shorten_loci,
                row_label_converter=loi.label_converter_shorten_loci,
                row_cluster=True, col_cluster=True, threads=threads,
                figsize=(key_side_size, key_side_size))

        # loci_key set:
        loci = []
        if chromhmm:
            loci.append(("chromhmm", 10))
        if default:
            loci.append(("default", 22))
        if repeats:
            loci.append(("repeats", 15))

        for loci_key, loci_side_size in loci:
            process_intersection_metric(
                loci_dict[key], loci_dict[loci_key],
                outdir / "{}@{}.csv".format(key, loci_key), pdf,
                adjustments=_adjustment() if loci_key not in {"default"} else _adjustment_wrc(),
                col_label_converter=loi.label_converter_shorten_loci,
                row_label_converter=loi.label_converter_shorten_loci,
                row_cluster=True, col_cluster=False, threads=threads,
                figsize=(loci_side_size, key_side_size))

            process_intersection_metric(
                loci_dict[loci_key], loci_dict[key],
                outdir / "{}@{}.csv".format(loci_key, key), pdf,
                adjustments=_adjustment() if key not in {"default"} else _adjustment_wrc(),
                col_label_converter=loi.label_converter_shorten_loci,
                row_label_converter=loi.label_converter_shorten_loci,
                row_cluster=False, col_cluster=True, threads=threads,
                figsize=(key_side_size, loci_side_size))

        # consensus
        consensus_paths = sorted(
            loci_dict['zinbra_{}'.format(consensus_type)] + loci_dict['golden_{}'.format(
                consensus_type)],
            key=lambda f: f.name
        )
        if consensus:
            process_intersection_metric(
                loci_dict[key], consensus_paths,
                outdir / "{}@{}.csv".format(key, consensus_type), pdf,
                adjustments=_adjustment(),
                col_label_converter=loi.label_converter_shorten_loci,
                row_label_converter=loi.label_converter_shorten_loci,
                row_cluster=True, col_cluster=False, threads=threads, figsize=(20, key_side_size))

        # YDS or ODS, but not "_ODS_without_YDS_median_consensus.bed"
        if consensus_yo:
            consensus_yo_paths\
                = [p for p in consensus_paths if "DS" in p.name and "without" not in p.name]
            process_intersection_metric(
                loci_dict[key], consensus_yo_paths,
                outdir / "{}@{}_yo.csv".format(key, consensus_type), pdf,
                adjustments=_adjustment(),
                col_label_converter=loi.label_converter_shorten_loci,
                row_label_converter=loi.label_converter_shorten_loci,
                row_cluster=True, col_cluster=False, threads=threads, figsize=(12, key_side_size))


def process_intersection_metric(a_paths, b_paths, df_path: Path, pdf, **kw):
    df = bm.load_or_build_metrics_table(a_paths, b_paths, df_path)

    # print to pdf:
    bm.plot_metric_heatmap("IM: {}".format(df_path.name), df, save_to=pdf, **kw)
    return df


def init_pdf_info(pdf):
    # TODO: titles

    d = pdf.infodict()
    d['Title'] = 'Report: Intersection metric at different loci'
    d['Author'] = 'JetBrains Research BioLabs'
    d['Subject'] = 'outliers'
    # d['Keywords'] = 'outliers jetbrains aging'
    d['CreationDate'] = datetime.datetime.today()
    d['ModDate'] = datetime.datetime.today()


def _pvalues_above_thr(thr005, thr001):
    def inner(label):
        anns = []
        if label in thr005:
            anns.append(("0.05", "blue"))
        else:
            anns.append(("0.05", "lightgray"))

        if label in thr001:
            anns.append(("0.01", "green"))
        else:
            anns.append(("0.01", "lightgray"))
        return anns

    return inner


def test_donors(tool, peaks_map, loci_dict, loci_key, outdir, threads, outliers_df):
    peaks_dict = peaks_map[tool]
    tool = tool or "all"  # presentable text for labels, we use 'None' for all loci

    loci_paths = loci_dict[loci_key]

    for hist in sorted(peaks_dict.keys()):
        peaks_paths = peaks_dict[hist]
        peaks_key = "{}_{}".format(tool, hist)

        result_plot_path = outdir / "plot_{}-{}_stat.pdf".format(peaks_key, loci_key)

        with PdfPages(str(result_plot_path)) as pdf:
            init_pdf_info(pdf)

            # Intersection metric: peaks@loci
            test_donors_by_metric(
                bm.load_or_build_metrics_table(peaks_paths, loci_paths,
                                               outdir / "{}@{}.csv".format(peaks_key, loci_key),
                                               threads=threads),
                hist, outliers_df, peaks_paths, pdf,
                outdir / "{}@{}_stat.csv".format(peaks_key, loci_key)
            )

            # Intersection metric: loci@peaks, e.g. for small loci, transpose to make plots
            # have donors at OY, loci at OX
            test_donors_by_metric(
                bm.load_or_build_metrics_table(loci_paths, peaks_paths,
                                               outdir / "{}@{}.csv".format(loci_key, peaks_key),
                                               threads=threads).T,
                hist, outliers_df, peaks_paths, pdf,
                outdir / "{}@{}_stat.csv".format(loci_key, peaks_key)
            )
    pass


def test_donors_by_metric(df, hist, outliers_df, peaks_paths, pdf, stats_df_path):
    ha = "two-sided"  # 'less', 'two-sided', or 'greater'

    ##########################################################################################
    # Stat test for each locus: Old vs Young
    if stats_df_path.exists():
        print("    Already exists, loading:", str(stats_df_path))
        loci_pvalues_df = pd.read_csv(stats_df_path, index_col=0)
    else:
        print("    Calculating:", str(stats_df_path))
        mask_od_group, mask_yd_group = split_by_age(hist, outliers_df, peaks_paths)
        df_ods = df[mask_od_group]
        df_yds = df[mask_yd_group]
        print("    Dfs: OD = {}, YD = {}".format(df_ods.shape, df_yds.shape))
        loci_pvalues_df = calc_loci_pvalues(df_ods, df_yds, ha)
        # Save results:
        loci_pvalues_df.to_csv(str(stats_df_path))

    # TODO: temp code: move to cacl section
    loci_pvalues_df = loci_pvalues_df.sort_values(by=["dfr_bh", "pvalue"])
    loci_pvalues_df.to_csv(str(stats_df_path))
    # print(loci_pvalues_df)

    # Plots

    # TODO: to one plot: pvalues + colors for adjusted threshold ?
    manhattan_plot(loci_pvalues_df.sort_values(by="pvalue"),
                   "pvalue",
                   "[{}] Mann whitney u test pvalues".format(stats_df_path.name),
                   correction="Uncorrected",
                   save_to=pdf,
                   xticks=len(loci_pvalues_df) < 100)

    manhattan_plot(loci_pvalues_df.sort_values(by="dfr_bh"),
                   "dfr_bh",
                   "[{}] Mann whitney u test pvalues".format(stats_df_path.name),
                   correction="Benjamini–Hochberg corrected",
                   save_to=pdf,
                   xticks=len(loci_pvalues_df) < 100)

    # Heatmap with selected cols:
    anns = [bm.color_annotator_age]
    if hist and outliers_df is not None:
        if hist in outliers_df.columns:
            anns.append(bm.color_annotator_outlier(outliers_df, hist))
    row_annotator = None if not anns else bm.color_annotator_chain(*anns)

    def loci_passed_thr(df, col, thr):
        return set(df.index[df[col] < thr].tolist())

    def plot_significant(col, title):
        thr001 = loci_passed_thr(loci_pvalues_df, col, 0.01)
        thr005 = loci_passed_thr(loci_pvalues_df, col, 0.05)
        thr01 = loci_passed_thr(loci_pvalues_df, col, 0.1)

        print("Loci passing", title, "threshold 0.1", col, len(thr01))
        print("Loci passing", title, "threshold 0.005", col, len(thr005))
        print("Loci passing", title, "threshold 0.001", col, len(thr001))
        if thr01:
            bm.plot_metric_heatmap(
                "IM {} with for {} < 0.1".format(stats_df_path.name, title),
                df.loc[:, thr01],
                save_to=pdf,
                adjustments=dict(left=0.15, top=0.95, right=0.9, bottom=0.3),
                row_cluster=False, col_cluster=False, figsize=(20, 15),
                col_color_annotator=_pvalues_above_thr(
                    {loi.label_converter_shorten_loci(s) for s in thr005},
                    {loi.label_converter_shorten_loci(s) for s in thr001}
                ),
                row_color_annotator=row_annotator,
                col_label_converter=loi.label_converter_shorten_loci,
                row_label_converter=bm.label_converter_donor_and_tool,
            )

    plot_significant("pvalue", "not-adjusted pvalues")
    plot_significant("dfr_bh", "adjusted pvalues (BH fdr)")


def calc_loci_pvalues(df_ods, df_yds, ha):
    pvalues = []

    columns = df_ods.columns.tolist()
    assert columns == df_yds.columns.tolist()
    n = len(columns)

    for i in range(n):
        try:
            pvalues.append(
                mannwhitneyu(df_ods.iloc[:, i], df_yds.iloc[:, i], alternative=ha).pvalue
            )
        except ValueError as e:
            print("Error: {} in file:\n{}".format(e, columns[i]))
            # TODO: problem_files.append(file)
            pvalues.append(np.nan)
    loci_pvalues_df = pd.DataFrame.from_dict(dict(loci=columns, pvalue=pvalues))
    loci_pvalues_df.index = loci_pvalues_df.loci
    loci_pvalues_df.drop("loci", inplace=True, axis=1)
    print("Not corrected pvalue, first 10 lowest pvalues:")
    print(loci_pvalues_df.sort_values(by="pvalue").head(10))
    # P-values correction
    #   see: http://www.statsmodels.org/dev/_modules/statsmodels/stats/multitest.html
    pvalues_col = loci_pvalues_df["pvalue"]
    pvalues_not_nan_mask = ~np.isnan(pvalues_col)
    _reject, pvals_corrected, *_ = multipletests(
        pvals=pvalues_col[pvalues_not_nan_mask],
        # fdr_bh, holm-sidak, bonferroni
        alpha=0.05, method="fdr_bh"
    )
    loci_pvalues_df["dfr_bh"] = np.nan
    loci_pvalues_df.loc[pvalues_not_nan_mask, "dfr_bh"] = pvals_corrected

    return loci_pvalues_df


def split_by_age(hist, outliers_df, peaks_paths):
    # Split: Old / Young donors
    donors_age_id = [loi.donor_order_id(p) for p in peaks_paths]
    groups = np.asarray([age for age, _id in donors_age_id])
    mask_od_group = groups == "OD"
    mask_yd_group = groups == "YD"
    print("    Peaks: [{}], OD: [{}], YD: [{}]".format(
        len(peaks_paths), np.sum(mask_od_group), np.sum(mask_yd_group)
    ))
    # Load Outliers info
    mask_not_outlier = None
    if hist and outliers_df is not None:
        if hist in outliers_df.columns:
            col = outliers_df.loc[:, hist]
            outliers_codes = [col["{}{}".format(age, id)] for age, id in donors_age_id]
            mask_not_outlier = np.asarray(outliers_codes) == 0
    if mask_not_outlier is None:
        print("    {}: No outliers info, use all donors".format(hist))
        mask_not_outlier = np.ones((len(peaks_paths), 1), dtype=bool)
    print("    Not outliers: [{}]".format(np.sum(mask_not_outlier)))
    # Young/old group without outliers
    mask_od_group = mask_od_group * mask_not_outlier
    mask_yd_group = mask_yd_group * mask_not_outlier
    print("    Peaks: OD w/o outliers: [{}], YD w/o outliers: [{}]".format(
        np.sum(mask_od_group),
        np.sum(mask_yd_group))
    )
    return mask_od_group, mask_yd_group


def manhattan_plot(pvalues_df, col_name, title, correction,
                   save_to=None, adjustments=None, xticks=True):
    plt.figure(figsize=(15, 8))
    n = pvalues_df.shape[0]

    ax = plt.subplot()
    ax.plot(range(n), 1/pvalues_df[col_name], marker=".", ls="")
    # ax.plot(range(n), 1/pvalues_df["pvalue"], marker=".", ls="")
    # ax.axhline(y=-np.log10(0.05), xmin=0, xmax=n, color="r", linestyle='dotted')
    ax.axhline(y=1/0.05, xmin=0, xmax=n, color="r", linestyle='dotted')
    ax.set_ylabel("{} pvalues (-log(p) scale )".format(correction))
    ax.set_yscale("log")
    ax.set_title("{} ({})".format(title, correction))
    if xticks:
        plt.xticks(range(n),
                   [loi.label_converter_shorten_loci(l) for l in pvalues_df.index],
                   rotation='vertical')

    adjustments = adjustments or {}
    plt.subplots_adjust(left=adjustments.get('left', 0.2),
                        right=adjustments.get('right', 0.8),
                        top=adjustments.get('top', 0.8),
                        bottom=adjustments.get('bottom', 0.3))
    bm.save_plot(save_to)


if __name__ == "__main__":
    parent_dir = os.path.dirname(os.path.realpath(__file__))
    project_root = os.path.abspath(os.path.join(parent_dir) + "/..")
    sys.path.insert(0, project_root)

    # Force matplotlib to not use any Xwindows backend.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    import reports.bed_metrics as bm
    import reports.loci_of_interest as loi

    from matplotlib.backends.backend_pdf import PdfPages

    _cli()
