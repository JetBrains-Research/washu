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


def bowtie_logs(folder):
    """Process bowtie logs processed by batch task"""
    # Here we rely on bowtie output
    READS = '.*reads processed: '
    REPORTED_ALIGNMENT = '.*reported alignment: '
    FAILED_TO_ALIGN = '.*failed to align: '
    SUPRESSED = '.*due to -m: '
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
            for line in open(dirpath + '/' + f, 'r'):
                if re.search(READS, line):
                    reads = re.sub(READS, '', line).strip()
                if re.search(REPORTED_ALIGNMENT, line):
                    aligned = re.sub(REPORTED_ALIGNMENT, '', line).strip()
                if re.search(FAILED_TO_ALIGN, line):
                    failed_to_align = re.sub(FAILED_TO_ALIGN, '', line).strip()
                if re.search(SUPRESSED, line):
                    supressed = re.sub(SUPRESSED, '', line).strip()
            df.loc[len(df)] = (f, reads, aligned, failed_to_align, supressed)
    return df


def process(folder):
    """Process bowtie logs and create summary report"""
    report = folder + '/bowtie_report.csv'
    bowtie_logs(folder).to_csv(report, index=False)
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
