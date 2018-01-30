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
import math
from collections import namedtuple

from downstream.aging import is_od, is_yd
from downstream.signals.signal_pca_fit_error_pvalue_permutation_test import collect_paths

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # nopep8

R2DistMetrics = namedtuple('R2DistMetrics', ['mean', 'median', 'p2', 'p5', 'p10', 'wdist'])
Record = namedtuple('Record', ['datatype', 'folder', 'norm', 'metrics'])


def _homogeneous_split(group_a, group_b):
    gr1 = set()
    for age_gr in [list(group_a), list(group_b)]:
        group_size = len(age_gr) // 2
        rnd.shuffle(age_gr)
        gr1.update(age_gr[0:group_size])

    return gr1


def _process(path: Path, simulations: int, seed: int, threads: int, plot=True) -> R2DistMetrics:
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
    half_ods = len(ods) // 2
    ods = ods[0:2 * half_ods]
    ods1 = ods[0:half_ods]
    ods2 = ods[half_ods:]

    yds = [c for c in df.columns if is_yd(c)]
    rnd.shuffle(yds)
    half_yds = len(yds) // 2
    yds = yds[0:2 * half_yds]
    yds1 = yds[0:half_yds]
    yds2 = yds[half_yds:]

    ##################################################
    # Signal DF:
    donors = sorted(ods + yds)
    signal = df.loc[:, donors]

    k = simulations // 3
    donors_groups_sizes = [k, k, k + (simulations % 3)]
    donors_groups_params = (ods, ods1, ods2), (yds, yds1, yds2), (donors, ods, yds)

    r2_list_path = path.with_suffix(".permutations.r2.csv")
    if not r2_list_path.exists():
        if threads == 1:
            r2_list = []
            for n, (group_ab, group_a, group_b) in zip(donors_groups_sizes, donors_groups_params):
                r2_list.extend(_multiple_homogeneous_split_r2(group_ab, group_a, group_b, signal,
                                                              n))
        else:
            with Pool(processes=threads) as pool:
                tasks = []
                for n, (group_ab, group_a, group_b) in zip(donors_groups_sizes,
                                                           donors_groups_params):
                    chunk_size = math.ceil(n / threads)
                    chunks = [min(n - start, chunk_size) for start in range(0, n, chunk_size)]

                    tasks.extend(
                        [pool.apply_async(_multiple_homogeneous_split_r2,
                                          (group_ab, group_a, group_b, signal, l)) for l in chunks]
                    )
                r2_list = list(chain(*(task.get(timeout=3600 * timeout_hours) for task in tasks)))

        # serialize:
        pd.DataFrame.from_dict({"r2": r2_list_path}).to_csv(
            str(r2_list_path),
            index=None
        )
    else:
        r2_list = pd.DataFrame.from_csv(str(r2_list_path), index_col=None).r2.tolist()
        assert len(r2_list) == simulations,\
            "Expected {} simulations, but was {} in: {}".format(simulations, len(r2_list),
                                                                r2_list_path)

    rr = np.asarray(r2_list)
    # Wasserstein distance, Earth mover's distance
    wdist = np.sqrt(np.mean((rr - 1) ** 2))

    dm = R2DistMetrics(np.mean(rr), np.median(rr), np.percentile(rr, 2), np.percentile(rr, 5),
                       np.percentile(rr, 10), wdist)

    print("mean = {}, median = {}, [min, max] = [{}, {}], [2%, 5%, 10%, 98%] = [{}, {}, {}, "
          "{}]".format(dm.mean, dm.median, np.min(rr), np.max(rr), dm.p2, dm.p5, dm.p10,
                       np.percentile(rr, 98)))

    if plot:
        plt.hist(rr, color="darkgray")
        plt.xlim(xmin=min(0, np.min(rr)), xmax=1)
        plt.title("R2 for {} hom-groups. Mean = {:.5f}, 50% = {:.5f}, 2% = {:.5f}".format(
            len(rr), dm.mean, dm.median, dm.p2)
        )
        plt.axvline(x=dm.mean, color="blue", label="R2 mean", linestyle="--", linewidth=0.9)
        plt.axvline(x=dm.median, color="black", label="R2 median", linestyle="--", linewidth=0.9)
        plt.axvline(x=dm.p2, color="red", label="R2 2% percentile", linestyle="--", linewidth=0.9)
        plt.xlabel("Each locus R2 for mean signal @ group1 vs group2")
        plt.legend()

        plt.savefig(str(path.with_suffix(".permutations.r2.png")))
        plt.close()

    return dm


def _homogeneous_split_r2(donors: List, group_a: List, group_b: List, signal: pd.DataFrame):
    gr1 = _homogeneous_split(group_a, group_b)
    gr1_means = signal.loc[:, gr1].mean(axis=1)
    gr2_means = signal.loc[:, list(set(donors) - gr1)].mean(axis=1)
    r2 = r2_score(gr1_means, gr2_means)
    return r2


def _multiple_homogeneous_split_r2(donors: List, group_a: List, group_b: List,
                                   signal: pd.DataFrame,
                                   simulations: int):
    return [_homogeneous_split_r2(donors, group_a, group_b, signal) for _i in range(simulations)]


def _cli():
    parser = argparse.ArgumentParser(
        description="For each given normalised signal file calculates r2 distribution for "
                    "donors random split in 2 groups of 3 kinds:"
                    "  * groups contain only old donors\n"
                    "  * groups contain only young donors\n"
                    "  * groups contain same numbers of old and young donors (i.e. mixed age).",
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
    process(paths, str(root / "report.permutation_r2.{}.csv".format(simulations)),
            seed, simulations, threads)


def process(paths: List[Path], output_path: str, seed: int, simulations: int, threads: int):
    n_paths = len(paths)
    records = []
    for i, path in enumerate(paths, 1):
        print("--- [{} / {}] -----------".format(i, n_paths))
        print("Process:", path)
        dm = _process(path, simulations=simulations, seed=seed, threads=threads)

        norm = re.sub('(.*_)|(\\.tsv$)', '', path.name)
        matches = re.match(".*/(H[a-z0-9]+)/.*", str(path), re.IGNORECASE)
        if matches:
            mod = matches.group(1)
        elif "/meth/" in str(path):
            mod = "meth"
        else:
            mod = "N/A"
        records.append(Record(mod, str(path.parent), norm, dm))

    # sort by first 3 cols: (mod, hist, norm)
    records.sort(key=lambda r: r[:3])
    if len(paths) > 1:
        print("====================")
        for i, r in enumerate(records, 1):
            print("{}. mean = {}, 50% = {}, 2% = {}, : [{}] {}".format(
                i, r.metrics.mean, r.metrics.mean, r.metrics.p2, r.norm, r.folder))

        df = pd.DataFrame.from_records(
            [(r.datatype, r.folder, r.norm, *r.metrics) for r in records],
            columns=["modification", "file", "normalization", *(R2DistMetrics._fields)]
        )
        # table: mod, folder, norm, error
        df.to_csv(output_path, index=None)

        print("R2 distribution features saved to:", output_path)


if __name__ == "__main__":
    _cli()
