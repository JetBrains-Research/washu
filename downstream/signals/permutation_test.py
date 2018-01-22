import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import numpy.random as rnd
from sklearn.decomposition import PCA
from multiprocessing import Pool

from downstream.aging import is_od, is_yd
from downstream.signals.signals_visualize import pca_separation_fit_error

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # nopep8


def _random_split_error(donors, ods, yds, x_r):
    gr1 = set()
    for age_gr in [list(ods), list(yds)]:
        group_size = len(age_gr) // 2
        rnd.shuffle(age_gr)
        gr1.update(age_gr[0:group_size])

    y = [0 if d in gr1 else 1 for d in donors]
    return pca_separation_fit_error(x_r, y)


def _process(path: Path, simulations: int, seed: int, threads: int, plot=True):
    # In case of 14 ODS, 14 YDS we expect about C[20,10]*C[20,10]/2 different
    # separation in random groups contained 10 ODS and 10 YDS donors
    # in each ~ 5.8 * 10^6 separations
    timeout_secs = 10
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

    pca = PCA(n_components=2)  # TODO?: as signals_visuzlize.py ?
    x_r = pca.fit_transform(signal.T.values)
    actual_error = pca_separation_fit_error(x_r,
                                            [0 if is_yd(d) else 1 for d in donors])

    if threads == 1:
        rnd_results = []
        for i in range(simulations):
            rnd_results.append(_random_split_error(donors, ods, yds, x_r))
    else:
        with Pool(processes=threads) as pool:
            tasks = [pool.apply_async(_random_split_error, (donors, ods, yds, x_r)) for _i in
                     range(simulations)]
            rnd_results = [task.get(timeout=timeout_secs) for task in tasks]

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
        plt.title("Pvalue for {} random age-independent groups = {}".format(len(rr), pvalue))
        plt.axvline(x=actual_error, color="red")
        plt.xlabel("PCA classification error")
        plt.savefig(str(path.with_suffix(".rnd_pvalue.png")))
        plt.close()

    return pvalue


def _cli():
    parser = argparse.ArgumentParser(
        description="For each given normalised signal file calculates pvalue for whether this "
                    "metric value could be achieved by chance in any age independent donors "
                    "split in groups. Each group contain same quantities of old/young donors.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("path",
                        help="A normalized signal *.tsv file or root folder containing such files")
    parser.add_argument('--seed', help="Random generator seed", type=int)
    parser.add_argument('-p', '--threads', help="Threads number for parallel processing",
                        type=int, default=4)
    parser.add_argument('-n', help="Simulations number to calculated pvalue", type=int,
                        default=100*1000)
    #
    args = parser.parse_args()
    root = Path(args.path)
    seed = args.seed
    simulations = args.n
    threads = args.threads

    print("Threads: {}, seed: {}, simulations: {}".format(threads, seed, simulations))

    if root.is_dir():
        paths = [str(p) for p in root.glob("**/*.filtered_*.tsv")]
        # TODO: hack, cleanup it
        paths = [p for p in paths if "/meth/" not in p]
    else:
        paths = [root]

    path2pvalue = []
    for path in paths:
        print("--------------")
        print("Process:", path)
        pvalue = _process(Path(path), simulations=simulations, seed=seed, threads=threads)
        path2pvalue.append((path, pvalue))
    path2pvalue.sort(key=lambda v: v[1])

    if len(paths) > 1:
        print("====================")
        for i, (path, pvalue) in enumerate(path2pvalue, 1):
            print("{}. {}: {}".format(i, pvalue, path))

        df = pd.DataFrame.from_records(path2pvalue, columns=["path", "pvalue"])
        out = str(root / "permutation.rnd_pvalue.csv")
        df.to_csv(out, index=None)

        print("P-values table saved to:", out)


if __name__ == "__main__":
    _cli()