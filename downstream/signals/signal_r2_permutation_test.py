import argparse
import re
from pathlib import Path
import pandas as pd
import numpy as np
import numpy.random as rnd
from multiprocessing import Pool
from sklearn.metrics import r2_score
from typing import List
from itertools import chain

from downstream.aging import is_od, is_yd
from downstream.signals.signal_pca_fit_error_pvalue_permutation_test import collect_paths

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # nopep8


def _homogeneous_split(ods, yds):
    gr1 = set()
    for age_gr in [list(ods), list(yds)]:
        group_size = len(age_gr) // 2
        rnd.shuffle(age_gr)
        gr1.update(age_gr[0:group_size])

    return gr1


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
    rnd.shuffle(ods)
    ods = ods[0:2 * (len(ods) // 2)]

    yds = [c for c in df.columns if is_yd(c)]
    rnd.shuffle(yds)
    yds = yds[0:2 * (len(yds) // 2)]

    ##################################################
    # Signal DF:
    donors = sorted(ods + yds)
    signal = df.loc[:, donors]

    if threads == 1:
        r2_list = _multiple_homogeneous_split_r2(donors, ods, yds, signal, simulations)
    else:
        chunk_size = simulations // threads
        chunks = [(i, min(simulations, i + chunk_size)) for i in range(0, simulations, chunk_size)]

        with Pool(processes=threads) as pool:
            tasks = [pool.apply_async(_multiple_homogeneous_split_r2,
                                      (donors, ods, yds, signal, e - s))
                     for s, e in chunks]
            r2_list = list(chain(*(task.get(timeout=3600*timeout_hours) for task in tasks)))

    rr = np.asarray(r2_list)
    r2_mean = np.mean(rr)
    r2_median = np.median(rr)

    print("mean = {}, median = {}, [min, max] = [{}, {}], [2%, 98%] = [{}, {}]".format(
        r2_mean, r2_median, np.min(rr), np.max(rr),
        np.percentile(rr, 2), np.percentile(rr, 98))
    )

    if plot:
        plt.hist(rr, color="darkgray")
        plt.xlim(xmin=0, xmax=1)
        plt.title("R2 for {} hom-groups. Mean = {:.5f}, Median = {:.5f}".format(
            len(rr), r2_mean, r2_median)
        )
        plt.axvline(x=r2_mean, color="red", label="R2 mean", linestyle="--", linewidth=0.9)
        plt.axvline(x=r2_median, color="black", label="R2 median", linestyle="--", linewidth=0.9)
        plt.xlabel("Each locus R2 for mean signal @ group1 vs group2")
        plt.legend()

        plt.savefig(str(path.with_suffix(".permutations.r2.png")))
        plt.close()

    return r2_mean, r2_median


def _homogeneous_split_r2(donors: List, ods: List, yds: List, signal: pd.DataFrame):
    gr1 = _homogeneous_split(ods, yds)
    gr1_means = signal.loc[:, gr1].mean(axis=1)
    gr2_means = signal.loc[:, list(set(donors) - gr1)].mean(axis=1)
    r2 = r2_score(gr1_means, gr2_means)
    return r2


def _multiple_homogeneous_split_r2(donors: List, ods: List, yds: List, signal: pd.DataFrame,
                                   simulations: int):
    return [_homogeneous_split_r2(donors, ods, yds, signal) for _i in range(simulations)]


def _cli():
    parser = argparse.ArgumentParser(
        description="For each given normalised signal file calculates r2 distribution for "
                    "donors random split in 2 groups, containing same fraction of old and "
                    "young donors.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("path",
                        help="A normalized signal *.tsv file or root folder containing such files")
    parser.add_argument('--seed', help="Random generator seed", type=int)
    parser.add_argument('-p', '--threads', help="Threads number for parallel processing",
                        type=int, default=4)
    parser.add_argument('-n', help="Simulations number to calculated pvalue", type=int,
                        default=100*1000)

    args = parser.parse_args()
    root = Path(args.path)
    seed = args.seed
    simulations = args.n
    threads = args.threads

    print("Threads: {}, seed: {}, simulations: {}".format(threads, seed, simulations))

    paths = collect_paths(root)
    process(paths, str(root / "report.permutation_r2.csv"), seed, simulations, threads)


def process(paths, output_path, seed, simulations, threads):
    n_paths = len(paths)
    records = []
    for i, path in enumerate(paths, 1):
        print("--- [{} / {}] -----------".format(i, n_paths))
        print("Process:", path)
        r2_mean, r2_median = _process(path, simulations=simulations, seed=seed, threads=threads)

        norm = re.sub('(.*_)|(\\.tsv$)', '', path.name)
        matches = re.match(".*/(H[a-z0-9]+)/.*", str(path), re.IGNORECASE)
        if matches:
            mod = matches.group(1)
        elif "/meth/" in str(path):
            mod = "meth"
        else:
            mod = "N/A"
        records.append((mod, str(path.parent), norm, r2_mean, r2_median))
    # sort by (median, mean)
    # records.sort(key=lambda v: v[-2:])
    # sort by (mod, hist, norm)
    records.sort(key=lambda v: (v[:3]))
    if len(paths) > 1:
        print("====================")
        for i, (_mod, path, norm, r2_mean, r2_median) in enumerate(records, 1):
            print("{}. mean = {}, median = {} : [{}] {}".format(i, r2_mean, r2_median, norm, path))

        df = pd.DataFrame.from_records(
            records,
            columns=["modification", "file", "normalization", "mean", "median"]
        )
        # table: mod, folder, norm, error
        df.to_csv(output_path, index=None)

        print("R2 distribution features saved to:", out)


if __name__ == "__main__":
    _cli()
