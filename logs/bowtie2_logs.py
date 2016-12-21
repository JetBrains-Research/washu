#!/usr/bin/env python
__author__ = 'oleg.shpynov@jetbrains.com'

import getopt
import sys
import pandas as pd
import os
import re

help_message = 'Script to process bowtie2 logs summary.'


def usage():
    print(help_message)


# Here we rely on bowtie2 output
BOWTIE2_READS = ' reads; of these:'
BOWTIE2_FAILED_TO_ALIGN = ' aligned 0 times'
BOWTIE2_ALIGNED = ' aligned exactly 1 time'
BOWTIE2_ALIGNED_MULTIPLE_TIMES = ' aligned >1 times'


def bowtie2_logs(folder):
    """Process bowtie2 logs processed by batch task"""
    print('Processing bowtie2 logs', folder)
    df = pd.DataFrame(columns=['sample', 'reads', 'aligned', 'not_aligned', 'supressed'])
    for dirpath, dirs, files in os.walk(folder):
        for f in files:
            if 'bowtie2' not in f or not re.search('.log$', f):
                continue
            reads = ''
            aligned = ''
            failed_to_align = ''
            supressed = ''
            with open(dirpath + '/' + f, 'r') as report:
                for line in report:
                    if re.search(BOWTIE2_READS, line):
                        reads = re.sub(BOWTIE2_READS, '', line).strip()
                    if re.search(BOWTIE2_ALIGNED, line):
                        aligned = re.sub(BOWTIE2_ALIGNED, '', line).strip()
                    if re.search(BOWTIE2_FAILED_TO_ALIGN, line):
                        failed_to_align = re.sub(BOWTIE2_FAILED_TO_ALIGN, '', line).strip()
                    if re.search(BOWTIE2_ALIGNED_MULTIPLE_TIMES, line):
                        supressed = re.sub(BOWTIE2_ALIGNED_MULTIPLE_TIMES, '', line).strip()
            df.loc[len(df)] = (f, reads, aligned, failed_to_align, supressed)
    return df


def process(folder):
    """Process bowtie logs and create summary report"""
    report = folder + '/bowtie2_report.csv'
    report_df = bowtie2_logs(folder)
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
