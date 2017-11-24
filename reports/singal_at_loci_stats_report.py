import os
import sys
import argparse
from collections import defaultdict
from itertools import chain

import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


def _cli():
    ########################################################################
    parser = argparse.ArgumentParser(
        description="Check whether loci changed or not by signal metrics (raw, rpm, "
                    "rpkm) for interesting loci sets",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--signal', metavar="PATH",
                        # /mnt/stripe/bio/experiments/signal_consensus_peaks/
                        # /mnt/stripe/bio/experiments/signal_loci_of_interest/
                        default="/mnt/stripe/bio/experiments/signal",
                        help="Processed signal dir")
    parser.add_argument('-o', '--out', required=True, metavar="PATH",
                        help="Output dir")
    args = parser.parse_args()
    signal_root = Path(args.signal)

    # /mnt/stripe/bio/experiments/aging/signal@loci_stats
    results_dir = Path(args.out)
    results_dir.mkdir(parents=True, exist_ok=True)
    ########################################################################

    # signal_dfs_by_datatype, signal_dfs_by_loci = build_signal_dfs(results_dir, signal_root)
    # print(signal_dfs_by_datatype[("H3K4me1", "rpkm")].head())

    stats_test(results_dir, signal_root)


def stats_test(results_dir, signal_root):
    signal_pvalues_df_path = results_dir / "signal_pvalues.csv"
    signal_adjusted_pvalues_df_path = results_dir / "signal_adjusted_pvalues.csv"

    if signal_pvalues_df_path.exists():
        signal_pvalues_df = pd.read_csv(signal_pvalues_df_path, index_col=0)
    else:
        signal_pvalues_df = calc_signal_pvalues(signal_pvalues_df_path, signal_root)

    print("Processed hypothesis:", len(signal_pvalues_df))
    print(signal_pvalues_df.head(10))
    signal_pvalues_df.index = signal_pvalues_df.name
    signal_pvalues_df.drop("name", inplace=True, axis=1)
    print("Not corrected pval, first 10 lowerest pvalues:")
    signal_pvalues_df["min"] = signal_pvalues_df.min(axis=1)
    print(signal_pvalues_df.sort_values(by="min").head(10).to_string(line_width=300))
    # P-values correction
    # see: http://www.statsmodels.org/dev/_modules/statsmodels/stats/multitest.html
    print("Adjust pvalues..")
    signal_pvalues_bh_df = signal_pvalues_df.copy().drop("min", axis=1)
    for c in signal_pvalues_bh_df.columns:
        pvalues = signal_pvalues_bh_df.loc[:, c]
        pvalues_not_nan_mask = ~np.isnan(pvalues)
        pvalues_not_nan = pvalues[pvalues_not_nan_mask]
        _reject, pvalues_corrected, *_ = multipletests(
            pvals=pvalues_not_nan,
            # fdr_bh, holm-sidak, bonferroni
            alpha=0.05, method="fdr_bh"
        )
        signal_pvalues_bh_df.loc[pvalues_not_nan_mask, c] = pvalues_corrected
    signal_pvalues_bh_df["min"] = signal_pvalues_bh_df.min(axis=1, skipna=True)
    signal_pvalues_bh_sorted_df = signal_pvalues_bh_df.sort_values(by="min")
    signal_pvalues_bh_sorted_df.to_csv(str(signal_adjusted_pvalues_df_path))

    # Passing FDR correction
    signal_pvalues_bh_sorted_df_005 = signal_pvalues_bh_sorted_df[
        signal_pvalues_bh_sorted_df["min"] < 0.05]
    print("Passing FDR 0.05 by any metric:", len(signal_pvalues_bh_sorted_df_005))
    # print(signal_pvalues_bh_sorted_df_005.head(10).to_string(line_width=300))
    print("Corrected, first 10 lowerest pvalues:")
    print(signal_pvalues_bh_sorted_df.head(10).to_string(line_width=300))
    print("Same records, but original pvalues:")
    print(signal_pvalues_df.loc[signal_pvalues_bh_sorted_df.head(10).index, :].to_string(
        line_width=300))

    # Plots:

    with PdfPages(str(results_dir / "signal_pvalues.pdf")) as pdf:
        for col in signal_pvalues_df.columns:
            loir.manhattan_plot(
                signal_pvalues_df.sort_values(by="min"), col,
                "Signal [{}] ODS vs YDS: Mann whitney u test pvalues".format(col),
                correction="Uncorrected",
                save_to=pdf
            )
            loir.manhattan_plot(
                signal_pvalues_bh_sorted_df.sort_values(by="min"), col,
                "Signal [{}] ODS vs YDS: Mann whitney u test pvalues".format(col),
                correction="Benjamini–Hochberg corrected",
                save_to=pdf
            )

            if (col != "min"):
                plot_signal_at_signif_loci("Uncorrected",
                                           signal_pvalues_df,
                                           col, pdf, signal_root)
                plot_signal_at_signif_loci("Benjamini–Hochberg corrected",
                                           signal_pvalues_bh_sorted_df, col, pdf, signal_root)


def plot_signal_at_signif_loci(title, df, col, pdf, signal_root):

    def loci_passed_thr(df, col, thr):
        return set(df.index[df[col] < thr].tolist())

    series_list = []
    thr01 = sorted(loci_passed_thr(df, col, 0.1))
    thr005 = sorted(loci_passed_thr(df, col, 0.05))
    thr001 = sorted(loci_passed_thr(df, col, 0.01))
    for label in thr01:
        datatype, loci = label.split('@')
        path = signal_root / datatype / loci / "{}_{}_data.csv".format(loci, col)
        series = pd.read_csv(path, index_col=0).iloc[:, 0]
        series.name = label
        series_list.append(series)

    if series_list:
        df = pd.DataFrame(series_list, ).T

        # Normalize by columns (by loci across all donors)
        df = ((df - np.min(df, axis=0)) / (np.max(df, axis=0) - np.min(df, axis=0)))

        row_annotator = bm.color_annotator_age

        bm.plot_metric_heatmap(
            "Average [{}] signal normalized in each loci to 0..1, significant {} pvalues < "
            "0.1".format(col, title),
            df,
            save_to=pdf,
            adjustments=dict(left=0.15, top=0.95, right=0.9, bottom=0.3),
            row_cluster=False, col_cluster=False, figsize=(20, 15),
            col_color_annotator=loir._pvalues_above_thr(
                {label_converter_shorten_loci(s) for s in thr005},
                {label_converter_shorten_loci(s) for s in thr001}
            ),
            row_color_annotator=row_annotator,
            col_label_converter=label_converter_shorten_loci,
            row_label_converter=bm.label_converter_donor_and_tool,
        )


def label_converter_shorten_loci(name):
    name = name.replace(".bed", "")
    name = name.replace("median_consensus", "mcs")
    name = name.replace("without", "w/o")
    return name


def calc_signal_pvalues(pvalues_df_path, signal_root):
    ha = "two-sided"  # 'less', 'two-sided', or 'greater'

    datatype_paths = [p for p in signal_root.iterdir() if p.is_dir()]

    loci_norm_paths_by_datatype, normalizations = collect_loci_normalizations(datatype_paths)
    print("Available normalizations", normalizations)

    signal_pvalues = defaultdict(list)
    missed_data = []
    problem_files = []
    for i, datatype_path in enumerate(datatype_paths, 1):
        data_type = datatype_path.name
        print("[{}/{}] Processing: {}".format(i, len(datatype_paths), data_type))

        loci_norm_paths = loci_norm_paths_by_datatype[data_type]

        # Collect missing data:
        for loci, dic in loci_norm_paths:
            missed_norms = normalizations - dic.keys()
            if missed_norms:
                missed_data.append("{}@{}:{}".format(loci, data_type, sorted(missed_norms)))

        # Test:
        for (loci, dic) in loci_norm_paths:
            pvalues = {}
            for norm, path in dic.items():
                df = pd.read_csv(path, index_col=0)
                df.index = [str.upper(s) for s in df.index]
                df.columns = [str.upper(s) for s in df.columns]
                df_ods = df.loc[[d for d in df.index if d.startswith("O")], :]
                df_yds = df.loc[[d for d in df.index if d.startswith("Y")], :]
                try:
                    pvalue = mannwhitneyu(df_ods.iloc[:, 0], df_yds.iloc[:, 0],
                                          alternative=ha).pvalue
                except ValueError as e:
                    print("Error: {} in file:\n{}".format(e, path))
                    problem_files.append(path)
                pvalues[norm] = pvalue

            signal_pvalues["name"].append("{}@{}".format(data_type, loci))
            for norm in normalizations:
                signal_pvalues[norm].append(pvalues.get(norm, np.nan))  # 1.0

    print("Missed files: ", len(missed_data))
    if missed_data:
        print("  first 10:", *missed_data[0:10])

    print("Errors occurred in files: ", len(problem_files))
    if problem_files:
        print("  first 10:", *problem_files[0:10])

    signal_pvalues_df = pd.DataFrame.from_dict(signal_pvalues)
    signal_pvalues_df.to_csv(str(pvalues_df_path))
    return signal_pvalues_df


def build_signal_dfs(signals_results_dir, signal_root):
    signal_dfs_by_datatype = {}
    signal_dfs_by_loci = {}

    series_by_loci_and_norm = defaultdict(lambda: defaultdict(list))

    datatype_paths = [p for p in signal_root.iterdir() if p.is_dir()]

    loci_norm_paths_by_datatype, normalizations = collect_loci_normalizations(datatype_paths)
    print("Available normalizations", normalizations)

    missed_data = []
    for i, datatype_path in enumerate(datatype_paths, 1):
        data_type = datatype_path.name
        print("[{}/{}] Processing: {}".format(i, len(datatype_paths), data_type))

        loci_norm_paths = loci_norm_paths_by_datatype[data_type]

        # Collect missing data:
        for loci, dic in loci_norm_paths:
            missed_norms = normalizations - dic.keys()
            if missed_norms:
                missed_data.append("{}@{}:{}".format(loci, data_type, sorted(missed_norms)))

        # To tables:
        datatype_series_by_norm = defaultdict(list)
        for (loci, dic) in loci_norm_paths:
            for norm, path in dic.items():
                series = pd.read_csv(path, index_col=0).iloc[:, 0]
                series.name = loci
                datatype_series_by_norm[norm].append(series)  # or (loci, series)

                series2 = series.copy()
                series2.name = data_type
                series_by_loci_and_norm[loci][norm].append(series2)

        for norm in normalizations:
            series_list = datatype_series_by_norm.get(norm, [])
            if series_list:
                df = pd.DataFrame(series_list, )
                df.to_csv(str(signals_results_dir / "signal_{}_{}.csv".format(data_type, norm)))
            else:
                df = None
            signal_dfs_by_datatype[(data_type, norm)] = df

    for loci, series_by_norm in series_by_loci_and_norm.items():
        for norm in normalizations:
            series = series_by_norm.get(norm, [])
            if series:
                df = pd.DataFrame(series, )
                df.to_csv(str(signals_results_dir / "signal_{}_{}".format(loci, norm)))
            else:
                df = None
            signal_dfs_by_loci[(loci, norm)] = df

    print("Missed files: ", len(missed_data))
    if (missed_data):
        print("  first 10:", *missed_data[0:10], sep="\n    ")
        print("Signal by datatype, missed records:",
              str(sum(1 for v in signal_dfs_by_datatype.values() if v is None)))
        print("  ", [k for k, v in signal_dfs_by_datatype.items() if v is None])
        print("Signal by loci, missed records:",
              str(sum(1 for v in signal_dfs_by_loci.values() if v is None)))
        # print("  ", [k for k, v in signal_dfs_by_loci.items() if v is None])

    return signal_dfs_by_datatype, signal_dfs_by_loci


def collect_loci_normalizations(datatype_paths):
    normalizations = {}
    loci_norm_paths_by_datatype = {}
    for i, datatype_path in enumerate(datatype_paths, 1):
        data_type = datatype_path.name
        print("[{}/{}] Collect paths: {}".format(i, len(datatype_paths), data_type))

        # Collect all loci x available normalizations info
        loci_norm_paths = []
        for loci_path in (p for p in datatype_path.iterdir() if p.is_dir()):
            loci = loci_path.name
            files = [p for p in loci_path.glob("**/*_data.csv")]

            available_norms = [
                f.name.replace("_data.csv", "").replace(loci + "_", "") for f in files
            ]
            loci_norm_paths.append((loci, dict(zip(available_norms, files))))
            # TODO cleanup:
            # break
        loci_norm_paths_by_datatype[data_type] = loci_norm_paths
        normalizations = set(chain(normalizations,
                                   *(dic.keys() for _loci, dic in loci_norm_paths)))
    return loci_norm_paths_by_datatype, normalizations


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
    import reports.loci_intersection_report as loir

    from matplotlib.backends.backend_pdf import PdfPages

    _cli()
