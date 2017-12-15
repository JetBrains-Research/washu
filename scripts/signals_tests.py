import numpy as np
import pandas as pd
import scipy
from scipy.stats import mannwhitneyu
from statsmodels.sandbox.stats.multicomp import multipletests

from scripts.util import *


def mann_whitney(x, y):
    try:
        return mannwhitneyu(x, y).pvalue
    except ValueError:
        return 1.0


def ttest(x, y):
    try:
        return scipy.stats.ttest_ind(x, y).pvalue
    except ValueError:
        return 1.0


def stat_test(f, test_name, test, fdr):
    print('Testing', test_name, f, 'fdr', fdr)
    df = pd.read_csv(f, sep='\t')
    # Drop contigs
    df = df.loc[[bool(re.match('chr[0-9XYM]+$', c)) for c in df['chr']]]
    ods = [c for c in df.columns.values if is_od(c)]
    yds = [c for c in df.columns.values if is_yd(c)]
    pvals = np.array([test(row[ods], row[yds]) for _, row in df.iterrows()])
    res = multipletests(pvals, fdr, "fdr_bh")
    h0_rejects = res[0]
    pvals_adj = res[1]
    df['pval'] = pvals
    df['pval_adj'] = pvals_adj
    df['od_mean'] = df[ods].mean(axis=1).to_frame('od_mean')['od_mean']
    df['yd_mean'] = df[yds].mean(axis=1).to_frame('yd_mean')['yd_mean']
    df['logfc'] = np.log(df['od_mean'] / df['yd_mean'])
    # Sort by adjusted pvalue
    pvals_adj_order = pvals_adj.argsort()
    df = df.loc[pvals_adj_order]
    h0_rejects = h0_rejects[pvals_adj_order]

    # Save results
    results = re.sub('\.tsv', '_{}.tsv'.format(test_name), f)
    df[['chr', 'start', 'end', 'yd_mean', 'od_mean', 'logfc', 'pval', 'pval_adj']] \
        .to_csv(results, sep='\t', index=None, header=True)
    print('Saved {} test results to {}'.format(sum(h0_rejects), results))

    # Save significant results
    if sum(h0_rejects) > 0:
        results_fdr = re.sub('\.tsv', '_{}_diff_fdr_{}.bed'.format(test_name, fdr), f)
        df.loc[h0_rejects][['chr', 'start', 'end']] \
            .to_csv(results_fdr, sep='\t', index=None, header=True)
        print('Saved {} significant results at FDR={} to {}'.format(sum(h0_rejects), fdr, results_fdr))


def process(work_dir, id):
    print('Stat testing', work_dir, id)
    id_dir = os.path.join(work_dir, id)
    for f in [f for f in os.listdir(id_dir) if
              re.match('.*{}_.*\\.tsv$'.format(id), f, flags=re.IGNORECASE)]:
        stat_test(os.path.join(id_dir, f), 'u', mann_whitney, 0.05)
        stat_test(os.path.join(id_dir, f), 't', ttest, 0.05)


def main():
    argv = sys.argv
    opts, args = getopt.getopt(argv[1:], "h", ["help"])
    if len(args) != 2:
        print("ARGUMENTS:  <work_dir> <id>\n"
              "CONVENTION: signal data is saved <folder>/<id>\n\n"
              "ARGS: " + ",".join(args))
        sys.exit(1)

    work_dir = args[0]
    id = args[1]

    process(work_dir, id)


if __name__ == "__main__":
    main()
