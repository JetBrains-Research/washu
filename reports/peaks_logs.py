#!/usr/bin/env python

import collections
import numpy as np
import getopt
import sys
import pandas as pd
import os
import re

__author__ = 'oleg.shpynov@jetbrains.com'
help_message = 'Script to process peaks RIP logs summary.'


def usage():
    print(help_message)


RipRecord = collections.namedtuple('RipRecord', ['file', 'peaks_file', 'reads', 'peaks', 'rip'])


def collect_rip_records(folder):
    """Collect all the Reads In Peaks records for given folder"""
    rips = []
    for dirpath, dirs, files in os.walk(folder):
        for f in files:
            if re.search('rip\\.csv$', f):
                with open(folder + '/' + f) as rip_file:
                    rips.append(RipRecord(*[line.rstrip('\n') for line in rip_file][1].split(',')))
    return rips


def report(folder):
    print('Process peaks logs processed by rip.sh script', folder)
    df = pd.DataFrame(np.empty((0,), dtype=[('sample', np.str),
                                            ('tags', np.int),
                                            ('peaks', np.int),
                                            ('rip', np.int),
                                            ('frip', np.int)]))
    rips = collect_rip_records(folder)
    for rr in rips:
        sample = rr.peaks_file.rpartition('/')[-1]
        reads = int(rr.reads)
        try:
            peaks = int(rr.peaks)
        except:
            peaks = 0
        try:
            rip = int(rr.rip)
        except:
            rip = 0
        frip = int(100 * rip / reads)
        df.loc[len(df)] = (sample, reads, peaks, rip, frip)
    return df.sort_values(by=['sample'])


def process_peaks_logs(folder, report_dir=None):
    """Process rip.csv files produced by rip.sh script and create summary"""
    df = report(folder)
    print(df)

    report_path = os.path.join(report_dir or folder, "peaks_report.csv")
    df.to_csv(report_path, index=False)
    print("Saved report", report_path)


def main():
    argv = sys.argv
    opts, args = getopt.getopt(argv[1:], "h", ["help"])
    if len(args) != 1:
        usage()
        sys.exit(1)
    process_peaks_logs(args[0])


if __name__ == "__main__":
    main()
