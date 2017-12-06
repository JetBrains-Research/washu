#!/usr/bin/env python
import getopt
import sys
import os
import re
import collections

import pandas as pd

from reports.peaks_logs import collect_rip_records, RipRecord

__author__ = 'oleg.shpynov@jetbrains.com roman.chernyatchik@jetbrains.com'

help_message = 'Script to process macs2 logs summary.'


def usage():
    print(help_message)


# Here we rely on macs2 output
MACS2_TAGS = '.*total tags in treatment:'
MACS2_REDUNDANT_RATE = '.*Redundant rate of treatment:'
MACS2_PAIRED_PEAKS = '.*paired peaks:'
MACS2_PREDICTED_FRAGMENT = '.*predicted fragment length is'
MACS2_ALTERNATIVE_FRAGMENTS = '.*alternative fragment length\(s\) may be'

MacsLogRecord = collections.namedtuple(
    'MacsLogRecord',
    ['sample', 'tags', 'redundant_rate', 'paired_peaks', 'fragment', 'alternatives']
)


def report(folder):
    print('Process macs2 logs processed by batch task', folder)

    macs_records = []
    for f in [f for f in os.listdir(folder)
              if re.match('.*macs2.*\\.log$', f, flags=re.IGNORECASE)]:
        tags = 0
        rr = 0.0
        paired_peaks = 0
        fragment = 0
        alt_fragments = ''
        with open(os.path.join(folder, f), 'r') as log:
            for line in log:
                if re.search(MACS2_TAGS, line):
                    tags = int(re.sub(MACS2_TAGS, '', line).strip())
                if re.search(MACS2_REDUNDANT_RATE, line):
                    rr = float(re.sub(MACS2_REDUNDANT_RATE, '',
                                      line).strip())
                if re.search(MACS2_PAIRED_PEAKS, line):
                    paired_peaks = int(re.sub(MACS2_PAIRED_PEAKS, '',
                                              line).strip())
                if re.search(MACS2_PREDICTED_FRAGMENT, line):
                    fragment = int(re.sub(MACS2_PREDICTED_FRAGMENT, '',
                                          line).replace('bps', '').strip())
                if re.search(MACS2_ALTERNATIVE_FRAGMENTS, line):
                    alt_fragments = re.sub(MACS2_ALTERNATIVE_FRAGMENTS, '',
                                           line).replace('bps', '').strip()

        macs_records.append(MacsLogRecord(f, tags, rr, paired_peaks,
                                          fragment, alt_fragments))

    df = pd.DataFrame(macs_records, columns=MacsLogRecord._fields)

    # Collect RIP records
    rips = collect_rip_records(folder)
    prefixes = (s.rpartition('_macs')[0] for s in df['sample'])
    sample_to_rips = [match_peaks_file(p, rips) for p in prefixes]

    df['peaks'] = [r.peaks for r in sample_to_rips]
    df['rip'] = [r.rip for r in sample_to_rips]
    df['frip'] = 100 * df['rip'] // df['tags']
    return df.sort_values(by=['sample'])


def match_peaks_file(x, records):
    """Find number of peaks in RipRecords"""
    filtered = [r for r in records if x in r.peaks_file]
    assert len(filtered) <= 1, "Ambiguous results for {}".format(x)
    return RipRecord("", "", 0, 0, 0) if (not filtered) else filtered[0]


def process_macs2_logs(folder, report_dir=None):
    """Process macs2 logs and create summary report"""
    df = report(folder)
    print(df)

    report_path = os.path.join(report_dir or folder, "macs2_report.csv")
    df.to_csv(report_path, index=False)
    print("Saved report", report_path)


def main():
    argv = sys.argv
    opts, args = getopt.getopt(argv[1:], "h", ["help"])
    if len(args) != 1:
        usage()
        sys.exit(1)
    process_macs2_logs(args[0])


if __name__ == "__main__":
    main()
