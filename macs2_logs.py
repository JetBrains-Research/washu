#!/usr/bin/env python
__author__ = 'oleg.shpynov@jetbrains.com'

import getopt
import sys
import pandas as pd
import os
import re

help_message = 'Script to process macs2 logs summary.'


def usage():
    print(help_message)


# Here we rely on macs2 output
MACS2_TAGS = '.*total tags in treatment:'
MACS2_REDUNDANT_RATE = '.*Redundant rate of treatment:'
MACS2_PAIRED_PEAKS = '.*paired peaks:'
MACS2_PREDICTED_FRAGMENT = '.*predicted fragment length is'
MACS2_ALTERNATIVE_FRAGMENTS = '.*alternative fragment length\(s\) may be'


def macs2_logs(folder):
    """Process macs2 logs processed by batch task"""
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
            with open(dirpath + '/' + f, 'r') as report:
                for line in report:
                    if re.search(MACS2_TAGS, line):
                        tags = int(re.sub(MACS2_TAGS, '', line).strip())
                    if re.search(MACS2_REDUNDANT_RATE, line):
                        rr = float(re.sub(MACS2_REDUNDANT_RATE, '', line).strip())
                    if re.search(MACS2_PAIRED_PEAKS, line):
                        paired_peaks = int(re.sub(MACS2_PAIRED_PEAKS, '', line).strip())
                    if re.search(MACS2_PREDICTED_FRAGMENT, line):
                        fragment = int(re.sub(MACS2_PREDICTED_FRAGMENT, '', line).replace('bps', '').strip())
                    if re.search(MACS2_ALTERNATIVE_FRAGMENTS, line):
                        alt_fragments = re.sub(MACS2_ALTERNATIVE_FRAGMENTS, '', line).replace('bps', '').strip()
            df.loc[len(df)] = (f, tags, rr, paired_peaks, fragment, alt_fragments)

    # Lines count data
    lcs = []
    # RiP data
    rips = []

    for dirpath, dirs, files in os.walk(folder):
        for f in files:
            # Peaks file
            if re.search('.(bed|broadPeak|narrowPeak)$', f):
                peaks = sum(1 for _ in open(folder + '/' + f))
                lcs.append((f, peaks))

            # _rip.txt file processing, see rip.sh
            if re.search('.txt$', f):
                with open(folder + '/' + f) as rip_file:
                    rips.append((f, int([line.rstrip('\n') for line in rip_file][0])))

    def lc_find(x):
        """Lines count"""
        rec = [lc[1] for lc in lcs if x.rpartition('_macs')[0] in lc[0]]
        return 0 if len(rec) == 0 else rec[0]

    def rip_find(x):
        """Read in Peaks"""
        rec = [rip[1] for rip in rips if x.rpartition('_macs')[0] in rip[0]]
        return 0 if len(rec) == 0 else rec[0]

    df['peaks'] = df['sample'].map(lambda x: lc_find(x))
    df['rip'] = df['sample'].map(lambda x: rip_find(x))
    df['frip'] = list(map(lambda x: int(x), 100 * df['rip'] / df['tags']))
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
