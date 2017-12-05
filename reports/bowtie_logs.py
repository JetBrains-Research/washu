#!/usr/bin/env python
import getopt
import sys
import pandas as pd
import os
import re

__author__ = 'oleg.shpynov@jetbrains.com'
help_message = 'Script to process bowtie logs summary.'


def usage():
    print(help_message)


# Here we rely on bowtie output
BOWTIE_READS = '.*reads processed: '
BOWTIE_REPORTED_ALIGNMENT = '.*reported alignment: '
BOWTIE_FAILED_TO_ALIGN = '.*failed to align: '
BOWTIE_SUPRESSED = '.*due to -m: '


def report(folder):
    """Process bowtie logs processed by batch task"""
    print('Processing bowtie logs', folder)
    df = pd.DataFrame(columns=['sample', 'reads', 'aligned', 'not_aligned', 'supressed'])
    for f in [f for f in os.listdir(folder) if re.match('.*bowtie.*\\.log$', f, flags=re.IGNORECASE)]:
        reads = ''
        aligned = ''
        failed_to_align = ''
        supressed = ''
        with open(os.path.join(folder, f), 'r') as log:
            for line in log:
                if re.search(BOWTIE_READS, line):
                    reads = re.sub(BOWTIE_READS, '', line).strip()
                if re.search(BOWTIE_REPORTED_ALIGNMENT, line):
                    aligned = re.sub(BOWTIE_REPORTED_ALIGNMENT, '',
                                     line).strip()
                if re.search(BOWTIE_FAILED_TO_ALIGN, line):
                    failed_to_align = re.sub(BOWTIE_FAILED_TO_ALIGN, '',
                                             line).strip()
                if re.search(BOWTIE_SUPRESSED, line):
                    supressed = re.sub(BOWTIE_SUPRESSED, '', line).strip()
        df.loc[len(df)] = (f, reads, aligned, failed_to_align, supressed)
    return df.sort_values(by=['sample'])


def process_bowtie_logs(folder, report_dir=None):
    """Process bowtie logs and create summary report"""
    df = report(folder)
    print(df)

    report_path = os.path.join(report_dir or folder, "bowtie_report.csv")
    df.to_csv(report_path, index=False)
    print("Saved report", report_path)


def main():
    argv = sys.argv
    opts, args = getopt.getopt(argv[1:], "h", ["help"])
    if len(args) != 1:
        usage()
        sys.exit(1)
    process_bowtie_logs(args[0])


if __name__ == "__main__":
    main()
