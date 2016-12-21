#!/usr/bin/env python
__author__ = 'oleg.shpynov@jetbrains.com'

import getopt
import sys
import pandas as pd
import os
import re

help_message = 'Script to process bowtie logs summary.'


def usage():
    print(help_message)


# Here we rely on bowtie output
BOWTIE_READS = '.*reads processed: '
BOWTIE_REPORTED_ALIGNMENT = '.*reported alignment: '
BOWTIE_FAILED_TO_ALIGN = '.*failed to align: '
BOWTIE_SUPRESSED = '.*due to -m: '


def bowtie_logs(folder):
    """Process bowtie logs processed by batch task"""
    print('Processing bowtie logs', folder)
    df = pd.DataFrame(columns=['sample', 'reads', 'aligned', 'not_aligned', 'supressed'])
    for dirpath, dirs, files in os.walk(folder):
        for f in files:
            if 'bowtie' not in f or not re.search('.log$', f):
                continue
            reads = ''
            aligned = ''
            failed_to_align = ''
            supressed = ''
            with open(dirpath + '/' + f, 'r') as report:
                for line in report:
                    if re.search(BOWTIE_READS, line):
                        reads = re.sub(BOWTIE_READS, '', line).strip()
                    if re.search(BOWTIE_REPORTED_ALIGNMENT, line):
                        aligned = re.sub(BOWTIE_REPORTED_ALIGNMENT, '', line).strip()
                    if re.search(BOWTIE_FAILED_TO_ALIGN, line):
                        failed_to_align = re.sub(BOWTIE_FAILED_TO_ALIGN, '', line).strip()
                    if re.search(BOWTIE_SUPRESSED, line):
                        supressed = re.sub(BOWTIE_SUPRESSED, '', line).strip()
            df.loc[len(df)] = (f, reads, aligned, failed_to_align, supressed)
    return df


def process(folder):
    """Process bowtie logs and create summary report"""
    report = folder + '/bowtie_report.csv'
    report_df = bowtie_logs(folder)
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
