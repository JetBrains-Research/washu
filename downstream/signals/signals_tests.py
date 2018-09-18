import numpy as np
import pandas as pd
import scipy
from scipy.stats import mannwhitneyu, wilcoxon
from statsmodels.sandbox.stats.multicomp import multipletests

from scripts.util import *
from downstream.aging import is_od, is_yd


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
    # Sort by pvalue
    pvals_order = pvals.argsort()
    df = df.loc[pvals_order]
    h0_rejects = h0_rejects[pvals_order]

    # Save results
    results = re.sub(r'\.tsv', '_{}.tsv'.format(test_name), f)
    df[['chr', 'start', 'end', 'yd_mean', 'od_mean', 'logfc', 'pval', 'pval_adj']] \
        .to_csv(results, sep='\t', index=None, header=True)
    print('Saved test results to', results)

    # Save significant results
    if sum(h0_rejects) > 0:
        results_fdr = re.sub(r'\.tsv', '_{}_diff_fdr_{}.bed'.format(test_name, fdr), f)
        df.loc[h0_rejects][['chr', 'start', 'end']] \
            .to_csv(results_fdr, sep='\t', index=None, header=True)
        print('Saved {} significant results at FDR={} to {}'.format(
            sum(h0_rejects), fdr, results_fdr))


def process(path):
    stat_test(path, 'u', mann_whitney, 0.05)
    stat_test(path, 'w', wilcoxon, 0.05)
    stat_test(path, 't', ttest, 0.05)


def main():
    argv = sys.argv
    opts, args = getopt.getopt(argv[1:], "h", ["help"])
    if len(args) != 1:
        print("ARGUMENTS:  <data.tsv>\n"
              "<data.tsv> - processed signal file to launch stat tests on"
              "ARGS: " + ",".join(args))
        sys.exit(1)

    path = args[0]
    process(path)


if __name__ == "__main__":
    main()
