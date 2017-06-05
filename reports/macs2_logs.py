#!/usr/bin/env python
import numpy as np

from reports.peaks_logs import collect_rip_records

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


def report(folder):
    print('Process macs2 logs processed by batch task', folder)
    df = pd.DataFrame(np.empty((0,), dtype=[('sample', np.str),
                                            ('tags', np.int),
                                            ('redundant_rate', np.float),
                                            ('paired_peaks', np.int),
                                            ('fragment', np.int),
                                            ('alternatives', np.str)]))
    for dirpath, dirs, files in os.walk(folder):
        for f in files:
            if 'macs' not in f or not re.search('.log$', f):
                continue
            tags = ''
            rr = ''
            paired_peaks = ''
            fragment = ''
            alt_fragments = ''
            with open(dirpath + '/' + f, 'r') as log:
                for line in log:
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

    # Collect RIP records
    rips = collect_rip_records(folder)
    df['peaks'] = df['sample'].map(lambda x: find_peaks(x.rpartition('_macs')[0], rips))
    df['rip'] = df['sample'].map(lambda x: find_rip(x.rpartition('_macs')[0], rips))
    df['frip'] = list(map(lambda x: int(x), 100 * df['rip'] / df['tags']))
    return df.sort_values(by=['sample'])


def find_peaks(x, rips):
    """Find number of peaks in RipRecords"""
    try:
        return int([rr.peaks for rr in rips if x in rr.peaks_file][0])
    except:
        return 0


def find_rip(x, rips):
    """Find Read in Peaks in RipRecords"""
    try:
        return int([rr.rip for rr in rips if x in rr.peaks_file][0])
    except:
        return 0


def process_macs2_logs(folder):
    """Process macs2 logs and create summary report"""
    path = folder + '/macs2_report.csv'
    df = report(folder)
    print(df)
    df.to_csv(path, index=False)
    print("Saved report", path)


def main():
    argv = sys.argv
    opts, args = getopt.getopt(argv[1:], "h", ["help"])
    if len(args) != 1:
        usage()
        sys.exit(1)
    process_macs2_logs(args[0])


if __name__ == "__main__":
    main()
