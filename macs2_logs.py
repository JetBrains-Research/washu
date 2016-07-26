#!/usr/bin/env python
__author__ = 'oleg.shpynov@jetbrains.com'

import getopt
import sys
import pandas as pd
import os
import re
import subprocess

help_message = 'Script to process macs2 logs summary.'


def usage():
    print(help_message)


def macs2_logs(folder):
    """Process macs2 logs processed by batch task"""
    # Here we rely on macs2 output
    TAGS = '.*total tags in treatment:'
    REDUNDANT_RATE = '.*Redundant rate of treatment:'
    PAIRED_PEAKS = '.*paired peaks:'
    PREDICTED_FRAGMENT = '.*predicted fragment length is'
    ALTERNATIVE_FRAGMENTS = '.*alternative fragment length\(s\) may be'
    print('Processing macs2 logs', folder)
    df = pd.DataFrame(columns=['sample', 'tags', 'redundant_rate', 'paired_peaks', 'fragment', 'alternatives'])
    for dirpath, dirs, files in os.walk(folder):
        for f in files:
            if 'macs' not in f or not re.search('.log$', f):
                continue
            tags = ''
            rr = ''
            paired_peaks = ''
            fragment = ''
            alt_fragments = ''
            for line in open(dirpath + '/' + f, 'r'):
                if re.search(TAGS, line):
                    tags = int(re.sub(TAGS, '', line).strip())
                if re.search(REDUNDANT_RATE, line):
                    rr = float(re.sub(REDUNDANT_RATE, '', line).strip())
                if re.search(PAIRED_PEAKS, line):
                    paired_peaks = int(re.sub(PAIRED_PEAKS, '', line).strip())
                if re.search(PREDICTED_FRAGMENT, line):
                    fragment = int(re.sub(PREDICTED_FRAGMENT, '', line).replace('bps', '').strip())
                if re.search(ALTERNATIVE_FRAGMENTS, line):
                    alt_fragments = re.sub(ALTERNATIVE_FRAGMENTS, '', line).replace('bps', '').strip()
            df.loc[len(df)] = (f, tags, rr, paired_peaks, fragment, alt_fragments)

    # Lines count data
    lcs = []
    # FRiP data
    rips = []
    for dirpath, dirs, files in os.walk(folder):
        for f in files:
            if re.search('.(bed|broadPeak|narrowPeak)$', f):
                lcs.append(
                    subprocess.Popen(['wc', '-l', folder + '/' + f], stdout=subprocess.PIPE).communicate()[0].decode('utf-8').strip())
            if re.search('.txt$', f):
                rips.append([line.rstrip('\n') for line in open(folder + '/' + f)][0])

    def lc_find(x):
        """Lines count"""
        rec = [lc.rpartition(' ')[0] for lc in lcs if x.rpartition('_macs')[0] in lc]
        return 0 if len(rec) == 0 else int(rec[0])

    def rip_find(x):
        """Read in Peaks"""
        rec = [rip.rpartition('\t')[2] for rip in rips if x.rpartition('_macs')[0] in rip]
        return 0 if len(rec) == 0 else int(rec[0])

    df['peaks'] = df['sample'].map(lambda x: lc_find(x))
    df['rip'] = df['sample'].map(lambda x: rip_find(x))
    df['frip'] = 100 * df['rip'] / df['tags']
    return df


def process(folder):
    """Process macs2 logs and create summary report"""
    report = folder + '/macs2_report.csv'
    report_df = macs2_logs(folder)
    print(report_df)
    report_df.to_csv(report, index=False)
    print("Saved report", report)


def main():
    argv = sys.argv
    opts, args = getopt.getopt(argv[1:], "h", ["help"])
    if len(args) != 1:
        usage()
        sys.exit(1)
    process(args[0])


if __name__ == "__main__":
    main()
