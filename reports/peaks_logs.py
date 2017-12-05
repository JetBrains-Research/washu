#!/usr/bin/env python

import collections
import getopt
import sys
import pandas as pd
import os
import re

__author__ = 'oleg.shpynov@jetbrains.com'
help_message = 'Script to process peaks RIP logs summary.'


def usage():
    print(help_message)


RipRecord = collections.namedtuple(
    'RipRecord', ['file', 'peaks_file', 'reads', 'peaks', 'rip']
)


def collect_rip_records(folder):
    """Collect all the Reads In Peaks records for given folder"""
    rips = []
    for f in [f for f in os.listdir(folder) if re.match('.*_rip\\.csv$', f, flags=re.IGNORECASE)]:
        with open(os.path.join(folder, f)) as rip_file:
            # skip header
            rip_file.readline()

            line = rip_file.readline().strip()
            records = line.split(',')
            assert len(records) == 5, \
                "Expected 5 comma separated values, but was {}: " \
                "line = '{}', file = {}".format(len(records), line, f)
            r = RipRecord(*records)

            # fix values if not set
            rips.append(r._replace(peaks=(int(r.peaks or 0)),
                                   rip=(int(r.rip or 0)),
                                   reads=int(r.reads)))
    return rips


# from collections import OrderedDict
def report(folder):
    print('Process peaks logs processed by rip.sh script', folder)

    records = collect_rip_records(folder)

    df = pd.DataFrame.from_items(
        [('sample', [rec.peaks_file.rpartition('/')[-1] for rec in records]),
         ('tags', [rec.reads for rec in records]),
         ('peaks', [rec.peaks for rec in records]),
         ('rip', [rec.rip for rec in records])]
    )
    df['frip'] = 100 * df['rip'] // df['tags']

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
