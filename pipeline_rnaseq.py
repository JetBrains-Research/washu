#!/usr/bin/env python
from pipeline_utils import *
import argparse

parser = argparse.ArgumentParser(description='RNA-seq data pipeline')
parser.add_argument(
    'path_to_directory', action=WritableDirectory, type=str,
    help='Path to directory with data to run pipeline')
parser.add_argument('path_to_indexes', action=WritableDirectory, type=str, help='Path to indexes')
parser.add_argument('genome', type=str, help='Genome')
args = parser.parse_args()

# Configuration
WORK_DIR = args.path_to_directory
GENOME = args.genome
INDEXES = os.path.join(args.path_to_indexes, GENOME)
CHROM_SIZES = os.path.join(INDEXES, GENOME + ".chrom.sizes")

print("WORK_DIR:", WORK_DIR)
print("GENOME:", GENOME)
print("INDEXES:", INDEXES)
print("CHROM_SIZES:", CHROM_SIZES)

##################
# Pipeline start #
##################

# Batch QC
run_bash("parallel/fastqc.sh", WORK_DIR)

# Batch STAR
run_bash("parallel/index_genome.sh", GENOME, INDEXES)
run_bash("parallel/index_star.sh", GENOME, INDEXES)
run_bash("parallel/star.sh", WORK_DIR, GENOME, os.path.join(INDEXES, 'star'))
move_forward(WORK_DIR, WORK_DIR + "_bams",
             ["*.bam", "*star*.log"])

WORK_DIR = WORK_DIR + "_bams"
os.chdir(WORK_DIR)

run_bash("parallel/bigwig.sh", CHROM_SIZES,
         os.path.join(INDEXES, 'star', GENOME + ".gtf"), WORK_DIR)
move_forward(WORK_DIR, WORK_DIR + "_bws",
             ["*.bw", "*bw.log", "*.bdg"])

# Batch RSEM
run_bash("parallel/index_rsem.sh", GENOME, INDEXES)
run_bash("parallel/rsem.sh", WORK_DIR, os.path.join(INDEXES, 'rsem', GENOME))
move_forward(WORK_DIR, WORK_DIR + "_rsem",
             ["*.results", "*.stat", "*rsem*", "genome_*.tsv", "transcriptome_*.tsv"])
