#!/usr/bin/env python
import os
import getopt
import re
import sys
import subprocess
import traceback

#################################################################
# Add project root folder to load path
from collections import namedtuple

this_file_path = os.path.realpath(__file__)
project_root_path \
    = os.path.abspath(os.path.join(os.path.dirname(this_file_path), os.pardir))
sys.path.insert(0, project_root_path)

from pipeline_utils import run_bash, move_forward  # nopep8
from reports.macs2_logs import process_macs2_logs  # nopep8

#################################################################

__author__ = 'oleg.shpynov@jetbrains.com'

help_message = '''
Usage:

python util.py find_input <file>
    Finds input given the file name. Heuristics: among all the files within\
     folder find file with "input" substring and
    most common subsequence with initial file.

python util.py macs_species <genome>
    Converts UCSC genome name to MACS.

python util.py effective_genome_fraction <genome> <chrom.sizes.path>
    Computes effective genome size, required for SICER.
'''


def usage():
    print(help_message)


#################################################################
# Age and donors utility code
#################################################################
Age = namedtuple('Age', 'name color prefix')
OLD = Age('O', 'blue', '')
YOUNG = Age('Y', 'red', '')


def is_od_input(c):
    return re.match('.*input.*od.*', str(c), flags=re.IGNORECASE) or \
           re.match('.*od.*input.*', str(c), flags=re.IGNORECASE)


def is_yd_input(c):
    return re.match('.*input.*yd.*', str(c), flags=re.IGNORECASE) or \
           re.match('.*yd.*input.*', str(c), flags=re.IGNORECASE)


def is_input(c):
    return is_od_input(c) or is_yd_input(c)


def is_od(c):
    return re.match('.*od\\d+.*', str(c), flags=re.IGNORECASE) and not is_input(c)


def is_yd(c):
    return re.match('.*yd\\d+.*', str(c), flags=re.IGNORECASE) and not is_input(c)


def is_od_or_yd(c):
    return re.match('.*[yo]d\\d+.*', str(c), flags=re.IGNORECASE) and not is_input(c)


def age(n):
    return re.search('[yo]d\\d+', str(n), flags=re.IGNORECASE).group(0)


def regions_extension(c):
    return re.match('.*(?:Peak|_peaks\.bed|island\.bed)$', str(c))
#################################################################


def lcs(x, y):
    """
    Finds longest common subsequence
    Code adopted from https://en.wikibooks.org/wiki/Algorithm_Implementation/
    Strings/Longest_common_subsequence#Python
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


def find_input(bam):
    filename = os.path.basename(bam)
    if 'input' in filename:
        return ''

    # Find all the files within folder
    dir_path = os.path.dirname(bam)
    f = []
    for (_, _, name) in os.walk(dir_path):
        f.extend(name)
        break

    def sort_function(x):
        return lcs(filename, x)

    inputs = [x for x in f if re.match('.*input.*\\.bam$', x)]
    if len(inputs) > 0:
        return max(inputs, key=sort_function)
    else:
        return ''


def macs_species(genome):
    'Convert genome to macs2 species encoding'
    if re.match('^hg[0-9]+$', genome):
        return 'hs'
    elif re.match('^mm[0-9]+$', genome):
        return 'mm'
    raise Exception('Unknown species {}'.format(genome))


def effective_genome_fraction(genome, chrom_sizes_path):
    """From MACS2 documentation:
    The default hs 2.7e9 is recommended for UCSC human hg18 assembly.
    Here are all precompiled parameters for effective genome size:
    hs: 2.7e9
    mm: 1.87e9
    ce: 9e7
    dm: 1.2e8"""
    chrom_length = int(run([['cat', chrom_sizes_path],
                            ['grep', '-v', 'chr_'],
                            ['awk', '{ L+=$2 } END { print L }']
                            ])[0].decode('utf-8').strip())
    if genome.startswith('mm'):
        size = 1.87e9
    elif genome.startswith('hg'):
        size = 2.7e9
    else:
        raise Exception('Unknown species {}'.format(genome))
    return size / chrom_length


def run_macs2(genome, chrom_sizes, name, *params, work_dirs):
    """
Defaults for MACS2 broad peak calling:
# effective genome size = 2.70e+09
# band width = 300
# model fold = [5, 50]
# qvalue cutoff for narrow/strong regions = 5.00e-02
# qvalue cutoff for broad/weak regions = 1.00e-01
# Larger dataset will be scaled towards smaller dataset.
# Range for calculating regional lambda is: 1000 bps and 10000 bps
# Broad region calling is on
# Paired-End mode is off

# Model params:
#  --bw, --mfold, --nolambda, --nomodel, --shift, --extsize
# Peak calling params:
#  --p, --q, --broad, --broad-cutoff
    """
    wd2result = {}
    skipped_wd = set()

    # Skip existing result folders: already calculated
    for wd in work_dirs:
        result_dir = '{}_macs2_{}'.format(wd, name)
        wd2result[wd] = result_dir

        if os.path.exists(result_dir):
            print('[Macs2] Already processed: ', result_dir)
            skipped_wd.add(wd)

    # Run MACS2
    unprocessed_workdirs = [wd for wd in work_dirs if wd not in skipped_wd]
    if unprocessed_workdirs:
        if os.path.exists(chrom_sizes):
            # -B produces bedgraph for signal
            params += ('-B',)

        run_bash("parallel/macs2.sh", genome, chrom_sizes, name,
                 "'{}'".format(" ".join([str(p) for p in params])),
                 *unprocessed_workdirs)

    # Move results to separate folder
    for wd in unprocessed_workdirs:
        result_dir = wd2result[wd]
        move_forward(
            wd, result_dir,
            ['*{}*'.format(name), '*rip.csv', '*_signal.bdg', '*_signal.bw'])

        process_macs2_logs(result_dir)

    return [os.path.join(wd, wd2result[wd]) for wd in work_dirs]


def run(commands, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    """Launches pipe of commands given stdin and final stdout, stderr"""

    # TODO: consider 'plumbum' library instead of this
    # https://plumbum.readthedocs.io/en/latest/local_commands.html

    processes = []
    _stdin = stdin

    for i, cmd in enumerate(commands):
        if i < len(commands) - 1:
            _stdout = subprocess.PIPE
            # Not clear how to collect stderr from chain, let's left
            # None here because result is more consistent:
            # * last cmd stderr is captured
            # * intermediate stderr is not missed, but not captured, goes to
            #    stderr
            # If you feel power, try to fix it + see tests
            _stderr = None
        else:
            _stdout = stdout
            _stderr = stderr

        p = subprocess.Popen(cmd, stdin=_stdin, stdout=_stdout,
                             stderr=_stderr)
        processes.append(p)
        _stdin = p.stdout

    for i in range(0, len(processes)):
        # noinspection PyBroadException
        try:
            if i < len(processes) - 1:
                # Allow p1 to receive a SIGPIPE if p2 exits.
                processes[i].stdout.close()
            else:
                sp = processes[i]
                out, err = sp.communicate()
                return out, err  # , sp.returncode

        except Exception:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            print("Error running: {}".format(commands))
            # exc_type below is ignored on 3.5 and later
            traceback.print_exception(exc_type, exc_value, exc_traceback,
                                      # limit=2,
                                      file=sys.stdout)


def main():
    argv = sys.argv
    opts, args = getopt.getopt(argv[1:], "h", ["help"])
    # Process help
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            return

    if len(args) == 2 and args[0] == 'find_input':
        print(find_input(args[1]))

    if len(args) == 2 and args[0] == 'macs_species':
        print(macs_species(args[1]))

    if len(args) == 3 and args[0] == 'effective_genome_fraction':
        print(effective_genome_fraction(args[1], args[2]))


if __name__ == "__main__":

    main()
