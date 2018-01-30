import argparse
import datetime
import sys
from collections import defaultdict
from itertools import chain
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


def _cli():
    data_root = Path("/mnt/stripe/bio")

    ########################################################################
    parser = argparse.ArgumentParser(
        description="Generates intersection reports for predefined loci sets and find loci"
                    "with significant intersection change between OLD and YOUNG donors group ("
                    "using Mann whitney u test)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-o', '--out', required=True, metavar="PATH",
                        help="Output dir")
    parser.add_argument('--loci', metavar="PATH",
                        # data_root / "experiments/loci_of_interest"
                        default=str(data_root / "raw-data/aging/loci_of_interest"),
                        help="Loci dir")
    parser.add_argument('-p', '--threads', type=int, default=4,
                        help="Threads number for parallel processing")
    parser.add_argument('--all', action="store_true",
                        help="Include failed tracks")
    parser.add_argument('--skip-plots', action="store_true",
                        help="Do not calc plots @ loci")
    parser.add_argument('--skip-stats', action="store_true",
                        help="Do not calc stats tests")
    parser.add_argument('--allpws', action="store_true",
                        help="Stat tests on all pathways")
    parser.add_argument('--tuned', action="store_true",
                        help="Use tuned peaks")
    parser.add_argument('--tool',
                        help="Process peaks only for these peak caller tool folder")
    parser.add_argument('--failed_tracks', metavar="PATH",
                        default="/mnt/stripe/bio/experiments/aging/Y20O20.failed_tracks.csv",
                        help="Failed tracks *.csv path")
    parser.add_argument('--peaks', metavar="PATH",
                        help="Custom peaks folder to use instead of predefined peaks list")
    args = parser.parse_args()

    threads = args.threads
    failed_tracks_df_path = args.failed_tracks
    exclude_failed_tracks = not args.all
    results_dir = Path(args.out)
    results_dir.mkdir(parents=True, exist_ok=True)
    loci_root = Path(args.loci)
    all_pathways = args.allpws
    skip_plots = args.skip_plots
    skip_stats = args.skip_stats
    custom_tool = args.tool
    ########################################################################

    loci_dict = loi.collect_loci(loci_root)

    if not args.peaks:
        if args.tuned:
            peaks_map = defaultdict(dict)
            bench_root = data_root / "experiments/configs/benchmark/benchmark"

            for hist_dir in (d for d in bench_root.glob("H*") if d.is_dir()):
                hist = hist_dir.name
                for tool_dir in (d for d in hist_dir.iterdir() if d.is_dir()):
                    tool = tool_dir.name

                    if custom_tool and tool != custom_tool:
                        continue

                    if tool == "macs_narrow":
                        # We decided to ignore Macs2 narrow
                        continue

                    peaks = loi.collect_peaks_in_folder(tool_dir)
                    if peaks:
                        peaks_map[tool][hist] = peaks

            tools_for_stat_test = sorted(peaks_map.keys())
        else:
            # TODO: legacy peaks, cleanup required
            # default peaks
            golden_peaks_root = data_root / "experiments/aging/peak_calling"
            zinbra_peaks_root = data_root / "experiments/configs/Y20O20{}/peaks".format(
                "" if exclude_failed_tracks else "_full"
            )
            # tuned peaks
            zinbra_peaks_root_tuned = data_root / "experiments/configs/Y20O20_full/benchmark_peaks"
            golden_peaks_root_tuned = data_root / "experiments/aging/peak_calling/benchmark_peaks"

            peaks_map = {
                "zinbra": loi._collect_zinbra_peaks(zinbra_peaks_root),
                "golden": loi._collect_golden_peaks(golden_peaks_root, exclude_failed_tracks),
                "zinbra_tuned": loi._collect_zinbra_peaks(zinbra_peaks_root_tuned),
                "golden_tuned": loi._collect_golden_peaks(golden_peaks_root_tuned, None)
            }
            tools_for_stat_test = ["zinbra_tuned", "golden_tuned"]
    else:
        tool = "tool"
        peaks_map = {tool: {
            "histmod": loi.collect_peaks_in_folder(Path(args.peaks))
        }}
        tools_for_stat_test = [tool]

    failed_tracks_df = None
    if failed_tracks_df_path:
        failed_tracks_df = pd.read_csv(failed_tracks_df_path, delimiter="\t", skiprows=1,
                                       index_col="donor")
    ########################################################################
    # NOTCH pathways as loci:
    notch_pathways = [p for p in (loci_root / "aging_pathways").glob('R-HSA-266082*.bed')]
    notch_pathways.extend(
        [p for p in (loci_root / "aging_pathways").glob('R-HSA-1912399*.bed')])
    notch_pathways.extend(
        [p for p in (loci_root / "aging_pathways").glob('R-HSA-264460*_cds.bed')])
    notch_pathways.sort(key=lambda p: p.name)
    loci_dict["notch_pathways"] = notch_pathways

    loci_dict["wo_pathways"] = sorted(set(
        chain(*[loci_dict[lt] for lt in loci_dict.keys() if lt and ("pathways" not in lt)])),
        key=lambda p: p.name
    )

    # Some predefined loci set:
    default_paths = []
    for key in ["top_level_paths", "enhancers", "regulatory", "repeats", 'chromhmm']:
        if key in loci_dict:
            default_paths.extend(loci_dict[key])
        else:
            print("Annotations not found:", str(loci_root / key), file=sys.stderr)
    loci_dict["default"] = sorted(default_paths, key=lambda p: p.name)

    plot_sizes = {
        "notch_pathways": 20,
        "aging_pathways": 200,
        "chromhmm": 10,
        "chromhmm2": 10,
        "genes": 10,
        "repeats": 15,
        "wo_pathways": 20,
    }

    for cons in [key for key in loci_dict if "consensus" in key]:
        loci_dict[cons + "_yo"] \
            = [p for p in loci_dict[cons] if "DS" in p.name and "without" not in p.name]
        plot_sizes[cons + "_yo"] = 10

        loci_dict[cons + "_common"] \
            = [p for p in loci_dict[cons] if "DS" not in p.name]
        plot_sizes[cons + "_common"] = 10

    # ########## For donors #############################################################
    if not skip_plots:
        for tool in sorted(peaks_map.keys()):
            for lt in ["default", "wo_pathways",
                       "median_consensus", "median_consensus_common",
                       "weak_consensus", "weak_consensus_common"]:

                if lt in loci_dict:
                    print("----- [Report]: Donors {}@{} ----".format(tool, lt))
                    report_donors(tool, peaks_map, loci_dict, lt, plot_sizes.get(lt, 20),
                                  results_dir, threads, failed_tracks_df)

        # ########## For loci #############################################################
        # If custom peaks folder, skip plots, calc only stat test
        if not args.peaks:
            loci = sorted({k for k in loci_dict if "pathways" not in k})
            for i, lt_a in enumerate(loci):
                for j, lt_b in enumerate(loci):
                    idx = i * len(loci) + j + 1
                    print("----- {}/{} [Report]: {}@{} ----".format(
                        idx, len(loci) * len(loci),
                        lt_a, lt_b)
                    )

                    report(lt_a, lt_b, loci_dict, results_dir, threads,
                           plot_sizes.get(lt_a, 20), plot_sizes.get(lt_b, 20))

            # Pathways:
            for lt_a in [key for key in ["aging_pathways", "notch_pathways"] if key in loci_dict]:
                for lt_b in [k for k in loci_dict if k and k.endswith("_consensus_yo")]:
                    print("----- [Report]: {}@{} ----".format(lt_a, lt_b))
                    report(lt_a, lt_b, loci_dict, results_dir, threads,
                           plot_sizes.get(lt_a, 20), plot_sizes.get(lt_b, 20))

    # ########## Stat tests #############################################################
    if not skip_stats:
        loci_sets = [
            ("wo_pathways",),
            ("wo_pathways", "aging_pathways")
        ]
        if all_pathways:
            # noinspection PyTypeChecker
            loci_sets.append(("wo_pathways", "aging_pathways", "other_pathways"))

        for i, loci_set in enumerate(loci_sets, 1):
            print("----- {}/{} [Stat tests]: {} ----".format(i, len(loci_sets), loci_set))

            for lt in [lt for lt in loci_set if lt not in loci_dict]:
                print("  [Skipped] No loci paths for: ", lt)
                continue

            test_donors(tools_for_stat_test, peaks_map, loci_set, loci_dict, results_dir,
                        failed_tracks_df, exclude_failed_tracks, threads)


def _adjustment():
    return dict(left=0.15, top=0.95, right=0.65, bottom=0.25)


def _adjustment_wrc():
    return dict(left=0.15, top=0.95, right=0.62, bottom=0.3)


def report_donors(tool, peaks_map, loci_dict, loci_key, key_side_size,
                  outdir, threads, failed_tracks_df):
    peaks_dict = peaks_map[tool]

    result_plot_path = outdir / "plot_{}@{}_by_donor.pdf".format(tool, loci_key)
    with PdfPages(str(result_plot_path)) as pdf:
        init_pdf_info(pdf)
        for hist in sorted(peaks_dict.keys()):
            anns = [bm.color_annotator_age]
            if hist and failed_tracks_df is not None:
                if hist in failed_tracks_df.columns:
                    anns.append(bm.color_annotator_outlier(failed_tracks_df, hist))
            annotator = None if not anns else bm.color_annotator_chain(*anns)

            process_intersection_metric(
                peaks_dict[hist], loci_dict[loci_key],
                outdir / "{}_{}@{}.csv".format(tool, hist, loci_key), pdf,
                adjustments=dict(left=0.15, top=0.95, right=0.9, bottom=0.3),
                row_cluster=False, col_cluster=False, threads=threads,
                figsize=(key_side_size, 15),
                row_color_annotator=annotator,
                col_label_converter=loi.label_converter_shorten_loci,
                row_label_converter=bm.label_converter_donor_and_tool,
            )


def report(a_key, b_key, loci_dict, outdir, threads, a_key_side, b_key_side):
    with PdfPages(str(outdir / "plot_{}@{}.pdf".format(a_key, b_key))) as pdf:
        init_pdf_info(pdf)

        process_intersection_metric(
            loci_dict[a_key], loci_dict[b_key],
            outdir / "{}@{}.csv".format(a_key, b_key),
            pdf,
            adjustments=_adjustment_wrc(),
            col_label_converter=loi.label_converter_shorten_loci,
            row_label_converter=loi.label_converter_shorten_loci,
            row_cluster=False, col_cluster=False, threads=threads,
            figsize=(a_key_side, b_key_side))

        process_intersection_metric(
            loci_dict[a_key], loci_dict[b_key],
            outdir / "{}@{}.csv".format(a_key, b_key),
            pdf,
            adjustments=_adjustment_wrc(),
            col_label_converter=loi.label_converter_shorten_loci,
            row_label_converter=loi.label_converter_shorten_loci,
            row_cluster=True, col_cluster=True, threads=threads,
            figsize=(a_key_side, b_key_side))


def process_intersection_metric(a_paths, b_paths, df_path: Path, save_to, threads=4, **kw):
    df = bm.load_or_build_metrics_table(a_paths, b_paths, df_path, threads=threads)

    # print to pdf:
    try:
        bm.plot_metric_heatmap("IM: {}".format(df_path.name), df, save_to=save_to, **kw)
    except ValueError as e:
        print(e, file=sys.stderr)

    return df


def init_pdf_info(pdf):
    # TODO: titles

    d = pdf.infodict()
    d['Title'] = 'Report: Intersection metric at different loci'
    d['Author'] = 'JetBrains Research BioLabs'
    d['Subject'] = 'failed_tracks'
    # d['Keywords'] = 'failed_tracks jetbrains aging'
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


def test_donors(tools, peaks_map, loci_set, loci_dict, outdir,
                failed_tracks_df, exclude_failed_tracks, threads):
    pvalue_dfs = []

    for lt in loci_set:
        loci_paths = loci_dict.get(lt, [])
        if not loci_paths:
            continue

        for tool in tools:
            peaks_dict = peaks_map[tool]

            for hist in sorted(peaks_dict.keys()):
                peaks_paths = peaks_dict[hist]
                peaks_key = "{}_{}".format(tool, hist)
                result_plot_path = outdir / "plot-stat_{}-{}.pdf".format(peaks_key, lt)

                with PdfPages(str(result_plot_path)) as pdf:
                    init_pdf_info(pdf)

                    if "pathways" not in lt:
                        # pathways are small, no sense in peaks@pathways
                        df = test_donors_by_metric(
                            bm.load_or_build_metrics_table(
                                peaks_paths, loci_paths,
                                outdir / "{}@{}.csv".format(peaks_key, lt),
                                threads=threads
                            ),
                            hist, failed_tracks_df, peaks_paths, pdf,
                            outdir / "stat-{}@{}.csv".format(peaks_key, lt),
                            exclude_failed_tracks,
                        )
                        df.index = ["{}@{}".format(peaks_key, s) for s in df.index]
                        pvalue_dfs.append(df)

                    # Intersection metric: loci@peaks, e.g. for small loci, transpose to make plots
                    # have donors at OY, loci at OX
                    df = test_donors_by_metric(
                        bm.load_or_build_metrics_table(
                            loci_paths, peaks_paths,
                            outdir / "{}@{}.csv".format(lt, peaks_key),
                            threads=threads
                        ).T,
                        hist, failed_tracks_df, peaks_paths, pdf,
                        outdir / "stat-{}@{}.csv".format(lt, peaks_key),
                        exclude_failed_tracks,
                    )
                    df.index = ["{}@{}".format(s, peaks_key) for s in df.index]
                    pvalue_dfs.append(df)

    # Summary report:
    loci_desc = "-".join(loci_set).replace("_pathways", "") + "_pathways"

    # sign, not_sign, all
    # TODO: remove duplicates
    loci_pvalues_df = pd.concat(pvalue_dfs)
    loci_pvalues_df.sort_values(by="pvalue", inplace=True)

    # P-values correction
    #   see: http://www.statsmodels.org/dev/_modules/statsmodels/stats/multitest.html
    pvalues_col = loci_pvalues_df["pvalue"]
    pvalues_not_nan_mask = ~np.isnan(pvalues_col)
    _reject, pvals_corrected, *_ = multipletests(
        pvals=pvalues_col[pvalues_not_nan_mask],
        # fdr_bh, holm-sidak, bonferroni
        alpha=0.05, method="fdr_bh"
    )
    loci_pvalues_df["fdr_bh"] = np.nan
    loci_pvalues_df.loc[pvalues_not_nan_mask, "fdr_bh"] = pvals_corrected
    loci_pvalues_df.sort_values(by=["fdr_bh", "pvalue"], inplace=True)
    all_pvalues_path = outdir / "stats.{}.all_pvalues.csv".format(loci_desc)
    loci_pvalues_df.to_csv(str(all_pvalues_path))
    print("Done. All pvalues: ", str(all_pvalues_path))

    pvalue001_df = loci_pvalues_df[loci_pvalues_df["pvalue"] < 0.01]
    pvalue001_df.to_csv(str(outdir / "stats.{}.notadjusted_pvalues0.01.csv".format(loci_desc)))
    print("Not adjusted pvalues < 0.01, first 10:")
    print(pvalue001_df.head(10).to_string(line_width=200, index=True))

    bh01_df = loci_pvalues_df[loci_pvalues_df["fdr_bh"] < 0.1]
    bh01_df.to_csv(str(outdir / "stats.{}.bh_pvalues0.1.csv".format(loci_desc)))
    print("BH adjusted pvalues, FDR < 0.1, first 10:")
    print(bh01_df.head(10).to_string(line_width=200, index=True))

    # # TODO: plots pvalues    for label in thr01:
    # def foo(df):
    #     for s in df.index:
    #         a, b = s.split('@', maxsplit=2)
    #         if a.endswith(".bed"):
    #             loci = a
    #         else:
    #             loci = b
    #         print()
    # # ['cd14_chromhmm18.hg19.17_ReprPCWk.bed', 'macs_broad_H3K27ac']
    # # ['macs_broad_H3K27ac', 'cd14_chromhmm18.hg19.17_ReprPCWk.bed']
    #
    #
    # #     datatype, loci = label.split('@')
    # #     path = signal_root / datatype / loci / "{}_{}_data.csv".format(loci, col)
    # #     series = pd.read_csv(path, index_col=0).iloc[:, 0]
    # #     series.name = label
    # #     series_list.append(series)


def test_donors_by_metric(df, hist, failed_tracks_df, peaks_paths, pdf, stats_df_path,
                          exclude_failed_tracks):
    ha = "two-sided"  # 'less', 'two-sided', or 'greater'

    ##########################################################################################
    # Stat test for each locus: Old vs Young
    if stats_df_path.exists():
        print("    Already exists, loading:", str(stats_df_path))
        loci_pvalues_df = pd.read_csv(stats_df_path, index_col=0)
    else:
        print("    Calculating:", str(stats_df_path))
        mask_od_group, mask_yd_group = split_by_age(hist, failed_tracks_df, peaks_paths,
                                                    exclude_failed_tracks)
        df_ods = df[mask_od_group]
        df_yds = df[mask_yd_group]
        print("    Dfs: OD = {}, YD = {}".format(df_ods.shape, df_yds.shape))
        loci_pvalues_df = calc_loci_pvalues(df_ods, df_yds, ha)
        # sort
        loci_pvalues_df.sort_values(by="pvalue", inplace=True)
        # Save results:
        print("    Done: ", str(stats_df_path))
        loci_pvalues_df.to_csv(str(stats_df_path))

    # Plots
    # TODO: to one plot: pvalues + colors for adjusted threshold ?
    manhattan_plot(loci_pvalues_df.sort_values(by="pvalue"),
                   "pvalue",
                   "[{}] Mann whitney u test pvalues".format(stats_df_path.name),
                   correction="Uncorrected",
                   save_to=pdf)

    _plot_donors_at_significant_loci(df, loci_pvalues_df, "pvalue", "not-adjusted pvalues",
                                     failed_tracks_df, hist, stats_df_path.name, pdf)

    return loci_pvalues_df


def _plot_donors_at_significant_loci(df,
                                     loci_pvalues_df, col,
                                     title, failed_tracks_df, hist, table_name,
                                     pdf):
    # Heatmap with selected cols:
    def loci_passed_thr(df, col, thr):
        return set(df.index[df[col] < thr].tolist())

    anns = [bm.color_annotator_age]
    if hist and failed_tracks_df is not None:
        if hist in failed_tracks_df.columns:
            anns.append(bm.color_annotator_outlier(failed_tracks_df, hist))
    row_annotator = None if not anns else bm.color_annotator_chain(*anns)

    thr001 = loci_passed_thr(loci_pvalues_df, col, 0.01)
    thr005 = loci_passed_thr(loci_pvalues_df, col, 0.05)
    thr01 = loci_passed_thr(loci_pvalues_df, col, 0.1)

    print("Loci passing", title, "threshold 0.1", col, len(thr01))
    print("Loci passing", title, "threshold 0.005", col, len(thr005))
    print("Loci passing", title, "threshold 0.001", col, len(thr001))
    if thr01:
        bm.plot_metric_heatmap(
            "IM {} with for {} < 0.1".format(table_name, title),
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

    return loci_pvalues_df


def split_by_age(hist, failed_tracks_df, peaks_paths, exclude_failed_tracks):
    # Split: Old / Young donors
    donors_age_id = [loi.donor_order_id(p) for p in peaks_paths]
    groups = np.asarray([age for age, _id in donors_age_id])
    mask_od_group = groups == "OD"
    mask_yd_group = groups == "YD"
    print("    Peaks: [{}], OD: [{}], YD: [{}]".format(
        len(peaks_paths), np.sum(mask_od_group), np.sum(mask_yd_group)
    ))
    # Load failed_tracks info
    mask_not_outlier = None
    if exclude_failed_tracks:
        if hist and failed_tracks_df is not None:
            if hist in failed_tracks_df.columns:
                col = failed_tracks_df.loc[:, hist]
                failed_tracks_codes = [col["{}{}".format(age, id)] for age, id in donors_age_id]
                mask_not_outlier = np.asarray(failed_tracks_codes) == 0

        if mask_not_outlier is None:
            print(
                "    {}: No failed tracks info, but exclude failed tracks option passed".format(
                    hist)
            )
            sys.exit(1)
    else:
        print("    {}: use all donors".format(hist))
        mask_not_outlier = np.ones((len(peaks_paths), 1), dtype=bool)

    print("    Not failed_tracks: [{}]".format(np.sum(mask_not_outlier)))
    # Young/old group without failed_tracks
    mask_od_group = mask_od_group * mask_not_outlier
    mask_yd_group = mask_yd_group * mask_not_outlier
    print("    Peaks: OD w/o failed_tracks: [{}], YD w/o failed_tracks: [{}]".format(
        np.sum(mask_od_group),
        np.sum(mask_yd_group))
    )
    return mask_od_group, mask_yd_group


def manhattan_plot(pvalues_df, col_name, title, correction,
                   save_to=None, adjustments=None,
                   xticks=None):

    if xticks is None:
        xticks = len(pvalues_df) < 100

    # TODO: --- hack!!! >>>>>>
    import matplotlib.pyplot as plt
    import downstream.loci_of_interest as loi
    import downstream.bed_metrics as bm
    # TODO: <<<<< hack!!! -----

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
    # Force matplotlib to not use any Xwindows backend.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    import downstream.bed_metrics as bm
    import downstream.loci_of_interest as loi

    from matplotlib.backends.backend_pdf import PdfPages

    import seaborn as sns
    sns.set(font_scale=0.7)
    _cli()
