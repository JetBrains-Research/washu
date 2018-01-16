#!/usr/bin/env python
import getopt
import sys
import pandas as pd
import os
import re

__author__ = 'oleg.shpynov@jetbrains.com'
help_message = 'Script to process bowtie2 logs summary.'


def usage():
    print(help_message)


# Here we rely on bowtie2 output
BOWTIE2_READS = ' reads; of these:'
BOWTIE2_FAILED_TO_ALIGN = ' aligned 0 times'
BOWTIE2_ALIGNED = ' aligned exactly 1 time'
BOWTIE2_ALIGNED_MULTIPLE_TIMES = ' aligned >1 times'


def report(folder):
    """Process bowtie2 logs processed by batch task"""
    print('Processing bowtie2 logs', folder)
    df = pd.DataFrame(columns=['sample', 'reads', 'aligned', 'not_aligned', 'supressed'])
    for f in [f for f in os.listdir(folder) if
              re.match('.*bowtie2.*\\.log$', f, flags=re.IGNORECASE)]:
        reads = ''
        aligned = ''
        failed_to_align = ''
        supressed = ''
        with open(os.path.join(folder, f), 'r') as log:
            for line in log:
                if re.search(BOWTIE2_READS, line):
                    reads = re.sub(BOWTIE2_READS, '', line).strip()
                if re.search(BOWTIE2_ALIGNED, line):
                    aligned = re.sub(BOWTIE2_ALIGNED, '', line).strip()
                if re.search(BOWTIE2_FAILED_TO_ALIGN, line):
                    failed_to_align \
                        = re.sub(BOWTIE2_FAILED_TO_ALIGN, '', line).strip()
                if re.search(BOWTIE2_ALIGNED_MULTIPLE_TIMES, line):
                    supressed = re.sub(BOWTIE2_ALIGNED_MULTIPLE_TIMES,
                                       '', line).strip()
        df.loc[len(df)] = (f, reads, aligned, failed_to_align, supressed)
    return df.sort_values(by=['sample'])


def process_bowtie2_logs(folder, report_dir=None):
    """Process bowtie logs and create summary report"""
    df = report(folder)
    print(df)

    report_path = os.path.join(report_dir or folder, "bowtie2_report.csv")
    df.to_csv(report_path, index=False)
    print("Saved report", report_path)


def main():
    argv = sys.argv
    opts, args = getopt.getopt(argv[1:], "h", ["help"])
    if len(args) != 1:
        usage()
        sys.exit(1)
    process_bowtie2_logs(args[0])


if __name__ == "__main__":
    main()
