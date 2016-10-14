#!/usr/bin/env python
__author__ = 'oleg.shpynov@jetbrains.com'

import getopt
import sys
import os

help_message = '''
    Script to find input given the file name.
    Heuristics: among all the files within folder find file with "input" substring and
    most common subsequence with initial file.
'''


def usage():
    print(help_message)


def lcs(x, y):
    """
    Finds longest common subsequence
    Code adopted from https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Longest_common_subsequence#Python
    """
    m = len(x)
    n = len(y)
    # An (m+1) times (n+1) matrix
    c = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if x[i - 1] == y[j - 1]:
                c[i][j] = c[i - 1][j - 1] + 1
            else:
                c[i][j] = max(c[i][j - 1], c[i - 1][j])

    def back_track(i, j):
        if i == 0 or j == 0:
            return ""
        elif x[i - 1] == y[j - 1]:
            return back_track(i - 1, j - 1) + x[i - 1]
        else:
            if c[i][j - 1] > c[i - 1][j]:
                return back_track(i, j - 1)
            else:
                return back_track(i - 1, j)

    return len(back_track(m, n))


def process(file):
    # Find all the files within folder
    dir_path = os.path.dirname(os.path.realpath(file))
    f = []
    for (_, _, name) in os.walk(dir_path):
        f.extend(name)
        break

    def sort_function(x):
        return lcs(file, x)

    extension = os.path.splitext(file)[1]
    inputs = [x for x in f if 'input' in x and extension in x]
    if len(inputs) > 0:
        print(max(inputs, key=sort_function))
    else:
        print('')


def main():
    # TODO[oleg] test me!
    if (lcs('UW_CD14_RO01746_k4me3_1_1_ENCFF001FYS.fastq', 'UW_CD14_input_ENCFF001HUV.fastq') <
            lcs('UW_CD14_RO01746_k4me3_1_1_ENCFF001FYS.fastq', 'Broad_CD14_2_input_ENCFF000CCW.fastq')):
        print("Regression found, add properly testing")
        sys.exit(1)

    argv = sys.argv
    opts, args = getopt.getopt(argv[1:], "h", ["help"])

    if len(args) != 1:
        usage()
        sys.exit(1)

    file = args[0]
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()

    process(file)


if __name__ == "__main__":
    main()