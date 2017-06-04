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
        if len(rr.peaks) > 0:
            peaks = int(rr.peaks)
        else:
            peaks = 0
        if len(rr.rip) > 0:
            rip = int(rr.rip)
        else:
            rip = 0
        frip = int(100 * rip / reads)
        df.loc[len(df)] = (sample, reads, peaks, rip, frip)
    return df.sort_values(by=['sample'])


def process_peaks_logs(folder):
    """Process rip.csv files produced by rip.sh script and create summary"""
    path = folder + '/peaks_report.csv'
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
    process_peaks_logs(args[0])


if __name__ == "__main__":
    main()
