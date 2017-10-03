#!/usr/bin/env python
from pipeline_utils import *
import argparse

parser = argparse.ArgumentParser(
    description='RNA-seq data pipeline for WashU cluster'
)
parser.add_argument(
    'path_to_directory', action=WritableDirectory, type=str,
    help='Path to directory with data to run pipeline')
parser.add_argument(
    'read_files', metavar='name read1 read2', type=str, nargs='+',
    help='reads files to be processed in format: name read1 read2 [name read1 '
         'read2 ...]\n'
         'read files must be absolute paths to read files')
args = parser.parse_args()

# Configuration
WORK_DIR = args.path_to_directory
GENOME = "hg19"  # Nothing to compare with aligned on hg38
INDEXES = "/scratch/artyomov_lab_aging/indexes/" + GENOME
CHROMSIZES = os.path.join(INDEXES, GENOME + ".chrom.sizes")

# Making soft links to raw data
os.chdir(WORK_DIR)
pairs = []
if len(args.read_files) % 3 != 0:
    raise argparse.ArgumentTypeError()
else:
    it = iter(args.read_files)
    triplets = zip(it, it, it)
    for name, read1, read2 in triplets:
        os.symlink(read1, name + "_1.fq.gz")
        os.symlink(read2, name + "_2.fq.gz")
        pairs.append(name + "_1.fq.gz")
        pairs.append(name + "_2.fq.gz")

# Batch QC
run_bash("parallel/fastqc.sh", WORK_DIR)

run_bash("parallel/star.sh", WORK_DIR, INDEXES, "\"" + " ".join(pairs) + "\"")
move_forward(WORK_DIR, WORK_DIR + "_bams",
             ["*.bam", "*star_*.log", "star_align_*"])
WORK_DIR = WORK_DIR + "_bams"
os.chdir(WORK_DIR)

run_bash("parallel/rnaseq_quality.sh", WORK_DIR)
move_forward(WORK_DIR, WORK_DIR + "_quality",
             ["*.rnastat", "rnaseq_quality.*"])

run_bash("parallel/bigwig.sh", CHROMSIZES, WORK_DIR)
move_forward(WORK_DIR, WORK_DIR + "_bws",
             ["*.bw", "*.bw.log", "*.bdg"])

run_bash("parallel/rsem.sh", WORK_DIR, INDEXES)
move_forward(WORK_DIR, WORK_DIR + "_rsem",
             ["*.results", "*.stat", "rsem_exp_*", "*.tsv"])
