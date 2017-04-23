#!/usr/bin/env python

"""
This is a chip-seq technical pipeline.

Usage:
  * Launch FastQC
  * Decide whether trimming or subsampling is required
  * Modify inplace copy of this pipeline
  * Launch pipeline, and wait for "Done" message

Conventions:
This pipeline uses folder naming as a steps, i.e. next step appends _suffix for the working folder,
stores results in new folder and change working folder if necessary.

Example:
run_6_7 -> run_6_7_trim -> run_6_7_trim_bams -> run_6_7_trim_bams_bws
                                             -> run_6_7_trim_bams_macs_0.01
                                             -> run_6_7_trim_bams_macs_broad_0.01

NOTE: python3 required
> source activate py3.5

author oleg.shpynov@jetbrains.com
"""
from reports.bowtie_logs import process_bowtie_logs
from pipeline_utils import *
from scripts.macs_util import run_macs2

parser = argparse.ArgumentParser(description='ULI ChIP-Seq data pipeline for WashU cluster')
parser.add_argument('path_to_directory', action=WritableDirectory, type=str,
                    help='Path to directory with data to run pipeline')
args = parser.parse_args()
#################
# Configuration #
#################
WORK_DIR = args.path_to_directory
GENOME = "hg19"
INDEXES = os.path.join("/scratch/artyomov_lab_aging/Y10OD10/chipseq/indexes", GENOME)
CHROM_SIZES = os.path.join(INDEXES, GENOME + ".chrom.sizes")
READS = 15  # Subsampling to 15mln reads

##################
# Pipeline start #
##################
print("Genomes and indices folder: ", INDEXES)
run_bash("index_genome.sh", GENOME, INDEXES)
run_bash("index_bowtie.sh", GENOME, INDEXES)

# Batch QC & multiqc
run_bash("fastqc.sh", WORK_DIR)
subprocess.run("multiqc " + WORK_DIR, shell=True)

# Batch Bowtie with trim 5 first base pairs
run_bash("index_bowtie.sh", GENOME, INDEXES)
run_bash("bowtie.sh", WORK_DIR, GENOME, INDEXES, "5")
WORK_DIR = move_forward(WORK_DIR, WORK_DIR + "_bams", ["*.bam", "*bowtie*.log"])
# multiqc is able to process Bowtie report
subprocess.run("multiqc " + WORK_DIR, shell=True)
# Create summary
process_bowtie_logs(WORK_DIR)

# Process insert size of BAM visualization
run_bash("fragments.sh", WORK_DIR)

# Batch BigWig visualization
run_bash("bigwig.sh", WORK_DIR, CHROM_SIZES)
move_forward(WORK_DIR, WORK_DIR + "_bws", ["*.bw", "*.bdg", "*bw.log"], copy_only=True)

# # Batch subsampling
# run_bash("subsample.sh", WORK_DIR, str(READS))
# WORK_DIR = move_forward(WORK_DIR, WORK_DIR + "_{}mln".format(READS), ["*{}*".format(READS)])
#
# # Batch BigWig visualization
# run_bash("bigwig.sh", WORK_DIR, CHROM_SIZES)
# move_forward(WORK_DIR, WORK_DIR + "_bws", ["*.bw", "*.bdg", "*bw.log"], copy_only=True)

########################
# Peak calling section #
########################

# Example for regular peak calling (https://github.com/taoliu/MACS)
run_macs2(WORK_DIR, GENOME, CHROM_SIZES, 'q0.01',
          '-q', 0.01)
# Example for broad peak calling (https://github.com/taoliu/MACS)
run_macs2(WORK_DIR, GENOME, CHROM_SIZES, 'broad_0.1',
          '--broad', '--broad-cutoff', 0.1)

# Default broad peak calling with modifications
run_macs2(WORK_DIR, GENOME, CHROM_SIZES, 'broad_0.1_mfold10-30_bw_300',
          '--broad', '--broad-cutoff', 0.1,
          '--mfold', 10, 30,
          '--bw', 300)
run_macs2(WORK_DIR, GENOME, CHROM_SIZES, 'broad_0.1_mfold10-30',
          '--broad', '--broad-cutoff', 0.1,
          '--mfold', 10, 30)
run_macs2(WORK_DIR, GENOME, CHROM_SIZES, 'broad_0.1_bw_300',
          '--broad', '--broad-cutoff', 0.1,
          '--bw', 300)
run_macs2(WORK_DIR, GENOME, CHROM_SIZES, 'broad_0.1_mfold2-100',
          '--broad', '--broad-cutoff', 0.1,
          '--mfold', 2, 100)
run_macs2(WORK_DIR, GENOME, CHROM_SIZES, 'broad_0.1_bw_150',
          '--broad', '--broad-cutoff', 0.1,
          '--bw', 150)
run_macs2(WORK_DIR, GENOME, CHROM_SIZES, 'broad_0.1_bw_600',
          '--broad', '--broad-cutoff', 0.1,
          '--bw', 600)


# Batch macs with different peak calling procedures settings
for P in [0.05]:
    run_macs2(WORK_DIR, GENOME, CHROM_SIZES, 'p{}'.format(P),
              '-p', P)
    run_macs2(WORK_DIR, GENOME, CHROM_SIZES, 'p{}_broad'.format(P),
              '-p', P, '--broad')
    run_macs2(WORK_DIR, GENOME, CHROM_SIZES, 'p{}_broad_nolambda'.format(P),
              '-p', P, '--broad', '--nolambda')
    # Cutoff for broad region.
    # If -p is set, this is a pvalue cutoff, otherwise, it's a qvalue cutoff. DEFAULT: 0.1
    CUTOFF = 0.5
    run_macs2(WORK_DIR, GENOME, CHROM_SIZES, 'p{}_broad_{}'.format(P, CUTOFF),
              '-p', P, '--broad', '--broad-cutoff', CUTOFF)

# Default is 0.01. For broad marks, you can try 0.05 as cutoff.
for Q in [0.01, 0.05]:
    run_macs2(WORK_DIR, GENOME, CHROM_SIZES, 'q{}'.format(Q),
              '-q', Q)
    run_macs2(WORK_DIR, GENOME, CHROM_SIZES, 'broad_{}'.format(Q),
              '--broad', '--broad-cutoff', Q)
    run_macs2(WORK_DIR, GENOME, CHROM_SIZES, 'broad_{}_nolambda'.format(Q),
              '--broad', '--broad-cutoff', Q, '--nolambda')
    run_macs2(WORK_DIR, GENOME, CHROM_SIZES, 'q{}_broad'.format(Q),
              '-q', Q, '--broad')
    # Cutoff for broad region.
    # If -p is set, this is a pvalue cutoff, otherwise, it's a qvalue cutoff. DEFAULT: 0.1
    CUTOFF = 0.5
    run_macs2(WORK_DIR, GENOME, CHROM_SIZES, 'q{}_broad_{}'.format(Q, CUTOFF),
              '-q', Q, '--broad', '--broad-cutoff', CUTOFF)

# Custom fragment peak calling option
# FRAGMENT = 138
# Q = 0.01
# run_macs2(WORK_DIR, GENOME, CHROM_SIZES, 'd{}_broad_{}'.format(FRAGMENT, Q),
#              '--broad', '--broad-cutoff', Q,
#              '--nomodel', '--shift', '0', '--extsize', FRAGMENT)

# # Batch macs14 with different peak calling procedures settings
# # P = 1e-5 is default for MACS14
# P = 0.00001
# NAME = '14_p{}'.format(P)
# FOLDER = '{}_macs_{}'.format(WORK_DIR, NAME)
# print(FOLDER)
# if not os.path.exists(FOLDER):
#     run_bash("macs14.sh", WORK_DIR, GENOME, str(P))
#     move_forward(WORK_DIR, FOLDER, ["*{}*".format(NAME)], copy_only=True)
#     process_macs2_logs(FOLDER)
#
# # Batch rseg
# rseg_suffix = '_rseg'
# if not os.path.exists(WORK_DIR + rseg_suffix):
#     run_bash("rseg.sh", WORK_DIR, GENOME, CHROM_SIZES)
#     move_forward(WORK_DIR, WORK_DIR + rseg_suffix, ["*-domains.bed", "*-scores.wig", "*-boundaries.bed",
#                                                     "*-boundary-scores.wig", "*-counts.bed"], copy_only=True)
