import argparse
import re
from multiprocessing import Pool
from typing import List
from itertools import chain
import math

from pathlib import Path
import pandas as pd
import numpy as np
import numpy.random as rnd
from statsmodels.stats.multitest import multipletests

from downstream.aging import is_od, is_yd
from downstream.signals.signals_visualize import pca_separation_fit_error, pca_signal

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # nopep8


def _random_split_error(donors, x_r):
    # Assume two groups: "group1" and "donors - group1"
    gr1 = set(rnd.choice(donors, len(donors) // 2, replace=False))
    y = [0 if d in gr1 else 1 for d in donors]
    return pca_separation_fit_error(x_r, y)


def _multiple_random_split_error(donors, x_r, simulations):
    return [_random_split_error(donors, x_r) for _i in range(simulations)]


def _process(path: Path, simulations: int, seed: int, threads: int, plot=True):
    # In case of 14 ODS, 14 YDS we expect about C[20,10]*C[20,10]/2 different
    # separation in random groups contained 10 ODS and 10 YDS donors
    # in each ~ 5.8 * 10^6 separations
    timeout_hours = 3
    rnd.seed(seed)  # is actual for tests in single thread mode

    ##################################################
    df = pd.read_csv(str(path), sep='\t')

    ##################################################
    # Make ODS, YDS groups even length for simplicity:
    ods = [c for c in df.columns if is_od(c)]
    yds = [c for c in df.columns if is_yd(c)]

    ##################################################
    # Signal DF:
    donors = sorted(ods + yds)
    signal = df.loc[:, donors]

    pca, x_r = pca_signal(signal)
    actual_error = pca_separation_fit_error(x_r,
                                            [0 if is_yd(d) else 1 for d in donors])

    if threads == 1:
        rnd_results = _multiple_random_split_error(donors, x_r, simulations)
    else:
        chunk_size = math.ceil(simulations / threads)
        chunks = [(i, min(simulations, i + chunk_size)) for i in range(0, simulations, chunk_size)]

        with Pool(processes=threads) as pool:
            tasks = [pool.apply_async(_multiple_random_split_error, (donors, x_r, e - s))
                     for s, e in chunks]
            rnd_results = list(chain(*(task.get(timeout=3600*timeout_hours) for task in tasks)))

    # hack not to get zero pvalue
    rnd_results.append(actual_error)
    rr = np.asarray(rnd_results)
    pvalue = np.sum(rr <= actual_error) / len(rr)

    print("[ACTUAL]: {}, [SIMUL]: [min, max] = [{}, {}], [2%, 98%] = [{}, {}]; 50% = {}, "
          "p-value: {}".format(actual_error,
                               np.min(rr[0:len(rr) - 1]), np.max(rr[0:len(rr) - 1]),
                               np.percentile(rr, 2), np.percentile(rr, 98), np.percentile(rr, 50),
                               pvalue
                               )
          )

    if plot:
        plt.hist(rr)
        plt.title("Pvalue for {} random groups = {}".format(len(rr), pvalue))
        plt.axvline(x=actual_error, color="red", label="ODS vs YDS error", linestyle="--",
                    linewidth=0.9)
        plt.xlabel("PCA classification error")
        plt.legend()

        plt.savefig(str(path.with_suffix(".permutations.pvalue.png")))
        plt.close()

    return pvalue, actual_error


def _cli():
    parser = argparse.ArgumentParser(
        description="For each given normalised signal file calculates pvalue for H0: pca "
                    "separation fit error for given ODS & YDS separation doesn't differ from  "
                    "random donors split into 2 groups.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("path",
                        help="A normalized signal *.tsv file or root folder containing such files")
    parser.add_argument('--seed', help="Random generator seed", type=int)
    parser.add_argument('-p', '--threads', help="Threads number for parallel processing",
                        type=int, default=4)
    parser.add_argument('-n', help="Simulations number to calculated pvalue", type=int,
                        default=100*1000)
    parser.add_argument('--fdr', help="Perform FDR control", action="store_true")
    parser.add_argument("--filter", required=True, metavar="SUBSTRINGS",
                        help="Comma separated file names filters.")

    args = parser.parse_args()
    root = Path(args.path)
    seed = args.seed
    simulations = args.n
    threads = args.threads
    fdr = args.fdr
    filters = [s.strip() for s in args.filter.split(",")]

    print("Threads: {}, seed: {}, simulations: {}".format(threads, seed, simulations))

    paths = [p for p in collect_paths(root) if all(f in p.name for f in filters)]

    if not paths:
        print("No suitable files found for '{}' and filters: {}".format(root, filters))
        exit(1)

    process(paths,
            str(root / "report.permutation_pvalue.{}.csv".format(simulations)),
            seed, simulations, threads, fdr)


def process(paths: List[Path], output_path: str, seed: int, simulations: int, threads: int,
            fdr: bool):
    n_paths = len(paths)
    records = []
    for i, path in enumerate(paths, 1):
        print("--- [{} / {}] -----------".format(i, n_paths))
        print("Process:", path)
        pvalue, actual_error = _process(path, simulations=simulations, seed=seed, threads=threads)

        norm = re.sub('(.*_)|(\\.tsv$)', '', path.name)
        matches = re.match(".*/(H[a-z0-9]+)/.*", str(path), re.IGNORECASE)
        if matches:
            mod = matches.group(1)
        elif "/meth/" in str(path):
            mod = "meth"
        else:
            mod = "N/A"
        records.append((mod, str(path.parent), norm, actual_error, pvalue))

    # sort by pvalue
    records.sort(key=lambda v: v[-1])
    if len(paths) > 1:
        print("====================")
        for i, (_mod, path, norm, actual_error, pvalue) in enumerate(records, 1):
            print("{}. {}: {} [{}] {}".format(i, pvalue, actual_error, norm, path))

        df = pd.DataFrame.from_records(
            records,
            columns=["modification", "file", "normalization", "error", "pvalue"]
        )

        if fdr:
            # FDR control:
            _reject, pvalues_corrected, *_ = multipletests(
                pvals=df["pvalue"],
                # fdr_bh, holm-sidak, bonferroni
                alpha=0.05, method="fdr_bh"
            )
            df["pvalue_corr"] = pvalues_corrected
            df.sort_values(by="pvalue_corr")

        # table: mod, folder, norm, error
        df.to_csv(output_path, index=None)

        print("P-values table saved to:", output_path)


def collect_paths(root: Path) -> List[Path]:
    if root.is_dir():
        # filter normalized signal files:
        paths = [p for p in root.glob("**/*.tsv") if p.name.startswith(p.parent.name + "_")]
        paths = [p for p in paths if not p.name.endswith(".bed.tsv")]
    else:
        paths = [root]
    return paths


if __name__ == "__main__":
    _cli()
