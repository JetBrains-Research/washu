__author__ = 'Roman Chernyatchik'
__email__ = 'roman.chernyatchik@jetbrains.com'

import argparse
from multiprocessing import Pool
from typing import List
from itertools import chain
import math
import datetime

from pathlib import Path
import pandas as pd
import numpy as np
import numpy.random as rnd
from statsmodels.stats.multitest import multipletests

from downstream.aging import is_od, is_yd
from downstream.signals.signals_util import extract_normalization, extract_datatype
from downstream.signals.signals_visualize import pca_separation_fit_error, pca_signal

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # nopep8


def _random_split_error(donors, x_r):
    # Assume two groups: "group1" and "donors - group1"
    gr1 = set(rnd.choice(donors, len(donors) // 2, replace=False))
    y = [0 if d in gr1 else 1 for d in donors]
    return pca_separation_fit_error(x_r, y)


def _multiple_random_split_error(donors, x_r, simulations, actual_error,
                                 opt_by_pvalue_cutoff, verbose):
    if not opt_by_pvalue_cutoff:
        return ([_random_split_error(donors, x_r) for _i in range(simulations)], False)
    else:
        res = []
        was_aborted = False
        # Check pvalue according threshold every 1000 simulations:
        chunk_size = 1000
        for fixed_size in [min(simulations - i, chunk_size) for i in range(0, simulations,
                                                                           chunk_size)]:
            res.extend(_random_split_error(donors, x_r) for _i in range(fixed_size))
            pvalue = _pvalue(actual_error, res, opt_by_pvalue_cutoff, simulations)
            if pvalue is None:
                # stop, no sense go further!
                was_aborted = True
                break
            elif verbose:
                print("###", len(res), pvalue, pvalue * len(res), ":", (simulations + 1) *
                      opt_by_pvalue_cutoff - 1)

        return (res, was_aborted)


def _process(path: Path, simulations: int, threads: int, opt_by_pvalue_cutoff=None,
             plot=True, verbose=True):

    # In case of 14 ODS, 14 YDS we expect about C[14,7]*C[14,7]/2 different
    # separation in random groups 'OY vs OY' contained 7 ODS and 7 YDS donors
    # Also C[14,7]/2 for "O vs O" and "Y vs Y"
    #
    # Total: (3432^2)/2 + 3432/2 + 3432/2 ~ 5.8 * 10^6 separations
    # Makes sense to reserve much less simulations number for "O vs O" and "Y vs Y".
    # E.g. min(1/3, 1500)

    df = pd.read_csv(str(path),
                     sep="," if path.suffix == ".csv" else "\t")

    if len(df) < 2:
        # cannot do PCA here
        print("[IGNORED] Not enough features {} < 2 in {}".format(len(df), path), flush=True)
        return path, 1, None

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

    rnd_results_path = path.with_suffix(".permutations.error.{}.csv".format(simulations))
    if not rnd_results_path.exists():
        if threads == 1:
            rnd_results, fake_pvalue = _multiple_random_split_error(
                donors, x_r, simulations, actual_error, opt_by_pvalue_cutoff, verbose
            )
            if verbose:
                print("[{}] fake pvalue: {}".format(path.name, fake_pvalue))
        else:
            fake_pvalue = False
            assert not opt_by_pvalue_cutoff,  \
                "Optimization by pvalue cutoff cannot be done in parallel for {}".format(path)

            timeout_hours = 3
            chunk_size = math.ceil(simulations / threads)
            chunks = [
                (i, min(simulations, i + chunk_size)) for i in range(0, simulations, chunk_size)
            ]

            with Pool(processes=threads) as pool:
                tasks = [pool.apply_async(_multiple_random_split_error,
                                          (donors, x_r, e - s, actual_error, None, verbose))
                         for s, e in chunks]
                rnd_results = list(
                    chain(*(task.get(timeout=3600*timeout_hours)[0] for task in tasks))
                )

        # serialize:
        if not fake_pvalue:
            pd.DataFrame.from_dict({"error": rnd_results}).to_csv(
                str(rnd_results_path),
                index=None
            )

    else:
        rnd_results = pd.read_csv(str(rnd_results_path), index_col=None).error.tolist()
        assert len(rnd_results) == simulations,\
            "Expected {} simulations, but was {} in: {}".format(simulations, len(rnd_results),
                                                                rnd_results_path)

    pvalue = _pvalue(actual_error, rnd_results, opt_by_pvalue_cutoff, simulations)
    if pvalue is None:
        pvalue = 1.0

    if verbose:
        rr = np.asarray(rnd_results)
        print("[ACTUAL]: {}, [SIMUL]: [min, max] = [{}, {}], [2%, 98%] = [{}, {}]; 50% = {}, "
              "p-value: {}".format(actual_error, np.min(rr), np.max(rr),
                                   np.percentile(rr, 2), np.percentile(rr, 98),
                                   np.percentile(rr, 50), pvalue)
              )

    if plot:
        plt.hist(np.asarray(rnd_results))
        plt.title("Pvalue for {} random groups = {}".format(len(rnd_results), pvalue))
        plt.axvline(x=actual_error, color="red", label="ODS vs YDS error", linestyle="--",
                    linewidth=0.9)
        plt.xlabel("PCA classification error")
        plt.legend()

        plt.savefig(str(path.with_suffix(".permutations.pvalue.{}.png".format(simulations))))
        plt.close()

    if not verbose:
        # print progress
        print('.', end="", flush=True)

    return path, pvalue, actual_error


def _pvalue(actual_error, errors, pvalue_cutoff, simulations):
    # hack not to get zero pvalue: pvalue = (better_errors number + 1) / (simulations number + 1)

    n_better_errors = sum(1 for e in errors if e <= actual_error)
    pvalue = (n_better_errors + 1) / (simulations + 1)

    if pvalue_cutoff and (pvalue > pvalue_cutoff):
        return None

    return pvalue


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

    parser.add_argument('--opt_max_pvalues',
                        help="Optimize permutation test so as to stop permutation test if pvalue "
                             "will not pass a threshold. By default p-values threshold is:"
                             "0.05 * opt_max_pvalues / hypothesis number. This threshold expects "
                             "to return <= opt_max_pvalues values passing FDR 0.05. All p-values "
                             "which won't pass the threshold are changed to 1.0. Could be useful "
                             "when default permutation test works too long. By default is turned "
                             "off.",
                        type=int,
                        default=None)
    parser.add_argument('--pvalue-thr', help="Custom p-value threshold", type=float)
    parser.add_argument('--fdr', help="Perform FDR control", action="store_true")
    parser.add_argument('--verbose', help="Detailed logging", action="store_true")
    parser.add_argument('--plots', help="Plot PCA classification error distribution for "
                                        "each loci",
                        action="store_true")
    parser.add_argument("--filter", metavar="SUBSTRINGS",
                        help="Comma separated file names filters.")

    args = parser.parse_args()
    root = Path(args.path)
    seed = args.seed
    simulations = args.n
    threads = args.threads
    fdr = args.fdr
    filters = [] if not args.filter else [s.strip() for s in args.filter.split(",")]
    opt_max_pvalues = args.opt_max_pvalues
    verbose = args.verbose
    plots = args.plots
    pvalue_threshold = args.pvalue_thr

    print("Threads: {}, seed: {}, simulations: {}".format(threads, seed, simulations))

    paths = [p for p in collect_paths(root) if all(f in p.name for f in filters)]

    if not paths:
        print("No suitable files found for '{}' and filters: {}".format(root, filters))
        exit(1)

    process(paths,
            str(root / "report.permutation_pvalue.{}.csv".format(simulations)),
            seed, simulations, threads, fdr, opt_max_pvalues=opt_max_pvalues,
            verbose=verbose, plot_simulations=plots,
            custom_pvalue_threshold=pvalue_threshold)


def process(paths: List[Path], output_path: str, seed: int, simulations: int, threads: int,
            fdr: bool, opt_max_pvalues=None, verbose=False, plot_simulations=True,
            custom_pvalue_threshold=None):

    started_ts = datetime.datetime.now()
    print("[{}] Starting..".format(started_ts.strftime('%Y-%m-%d %H:%M:%S')))

    if seed is not None:
        rnd.seed(seed)  # is actual for tests in single thread mode

    n_paths = len(paths)
    results = []

    if opt_max_pvalues is not None:
        timeout_hours = 10

        # According to BH fdr correction we need largest k: pk <= alpha * k / N
        # where pk is k-th of sorted pvalues, N - number of pvalues
        # So if we fix 'k' as max desired results number, we can estimate pvalue cutoff.
        # If we consider stricter FDR this cutoff wan't affect results. So let's take
        # rather weak FDR, e.g. 0.05 or 0.1
        pval_cutoff = custom_pvalue_threshold if custom_pvalue_threshold \
            else 0.05 * opt_max_pvalues / len(paths)
        print("P-value cutoff: ", pval_cutoff)

        # parallelize by loci
        with Pool(processes=threads) as pool:
            tasks = [pool.apply_async(_process,
                                      (p, simulations, 1, pval_cutoff),
                                      dict(plot=plot_simulations, verbose=verbose))
                     for p in paths]
            results = list(task.get(timeout=3600*timeout_hours) for task in tasks)
    else:
        # parallelize by simulations
        progress_step = n_paths // (100/5)  # every 5%
        for i, path in enumerate(paths, 1):
            if verbose:
                print("[{} / {}] Process: {}".format(i, n_paths, path))

            results.append(_process(
                path, simulations, threads, opt_by_pvalue_cutoff=None, plot=plot_simulations,
                verbose=verbose
            ))

            if not verbose and (i % progress_step == 0):
                print("  {}%", i // progress_step)

    records = []

    for (path, pvalue, actual_error) in results:
        norm = extract_normalization(path)
        dtype = extract_datatype(path)
        records.append((dtype, str(path.parent), norm, actual_error, pvalue))

    finished_ts = datetime.datetime.now()
    print("[{}] Pvalues calculated in {} s".format(finished_ts.strftime('%Y-%m-%d %H:%M:%S'),
                                                   finished_ts - started_ts))

    # sort by pvalue
    records.sort(key=lambda v: v[-1])
    if len(paths) > 1:
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
