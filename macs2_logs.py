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
                    tags = re.sub(TAGS, '', line).strip()
                if re.search(REDUNDANT_RATE, line):
                    rr = re.sub(REDUNDANT_RATE, '', line).strip()
                if re.search(PAIRED_PEAKS, line):
                    paired_peaks = re.sub(PAIRED_PEAKS, '', line).strip()
                if re.search(PREDICTED_FRAGMENT, line):
                    fragment = re.sub(PREDICTED_FRAGMENT, '', line).replace('bps', '').strip()
                if re.search(ALTERNATIVE_FRAGMENTS, line):
                    alt_fragments = re.sub(ALTERNATIVE_FRAGMENTS, '', line).replace('bps', '').strip()
            df.loc[len(df)] = (f, tags, rr, paired_peaks, fragment, alt_fragments)
    # Collect line sizes of peaks file
    wcs = []
    for dirpath, dirs, files in os.walk(folder):
        for f in files:
            if not re.search('.(bed|broadPeak|narrowPeak)$', f):
                continue
            wcs.append(subprocess.Popen(['wc', '-l', folder + '/' + f], stdout=subprocess.PIPE).communicate()[0].decode(
                'utf-8').strip())
    df['peaks'] = df['sample'].map(lambda x: [wc.rpartition(' ')[0] for wc in wcs if x.rpartition('_macs')[0] in wc][0])
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
