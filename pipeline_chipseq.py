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
from reports.macs2_logs import process_macs2_logs
from reports.bowtie_logs import process_bowtie_logs
from pipeline_utils import *

parser = argparse.ArgumentParser(description='ULI ChIP-Seq data pipeline for WashU cluster')
parser.add_argument('path_to_directory', action=WritableDirectory, type=str,
                    help='Path to directory with data to run pipeline')
args = parser.parse_args()

# Configuration
WORK_DIR = args.path_to_directory
GENOME = "hg19"
INDEXES = os.path.join("/scratch/artyomov_lab_aging/Y10OD10/chipseq/indexes", GENOME)
CHROM_SIZES = os.path.join(INDEXES, GENOME + ".chrom.sizes")
# READS = 15  # Subsampling to 15mln reads

print("Genomes and indices folder: ", INDEXES)
run_bash("index_genome.sh", GENOME, INDEXES)
run_bash("index_bowtie.sh", GENOME, INDEXES)

# Batch QC & multiqc
run_bash("fastqc.sh", WORK_DIR)
subprocess.run("multiqc " + WORK_DIR, shell=True)

# Batch Bowtie with trim 5 first base pairs
run_bash("index_bowtie.sh", GENOME, INDEXES)
run_bash("bowtie.sh", WORK_DIR, GENOME, INDEXES, "5")
WORK_DIR = move_forward(WORK_DIR, "_bams", ["*.bam", "*bowtie*.log"])
# multiqc is able to process Bowtie report
subprocess.run("multiqc " + WORK_DIR, shell=True)
# Create summary
process_bowtie_logs(WORK_DIR)

# Process insert size of BAM visualization
run_bash("fragments.sh", WORK_DIR)

# Batch BigWig visualization
run_bash("bigwig.sh", WORK_DIR, CHROM_SIZES)
move_forward(WORK_DIR, "_bws", ["*.bw", "*.bdg", "*bw.log"], copy_only=True)

# # Batch subsampling
# run_bash("subsample.sh", WORK_DIR, str(READS))
# WORK_DIR = move_forward(WORK_DIR, "_{}mln".format(READS), ["*{}*".format(READS)])
#
# # Batch BigWig visualization
# run_bash("bigwig.sh", WORK_DIR, CHROM_SIZES)
# move_forward(WORK_DIR, "_bws", ["*.bw", "*.bdg", "*bw.log"], copy_only=True)

# Batch macs with different peak calling procedures settings
P = 0.05
# macs_p_suffix = "_macs_p{}".format(P)
# if not os.path.exists(WORK_DIR + macs_p_suffix):
#     NAME = 'p{}'.format(P)
#     run_bash("macs2.sh", WORK_DIR, GENOME, CHROM_SIZES, NAME, '-B', '-p', str(P))
#     move_forward(WORK_DIR, macs_p_suffix, ["*{}*".format(NAME), '*.bw', '*.bdg'], copy_only=True)
#     process_macs2_logs(WORK_DIR + macs_p_suffix)

macs_broad_p_suffix = "_macs_broad_p{}".format(P)
if not os.path.exists(WORK_DIR + macs_broad_p_suffix):
    NAME = 'broad_p{}'.format(P)
    run_bash("macs2.sh", WORK_DIR, GENOME, CHROM_SIZES, NAME, '-B', '-p', str(P), '--broad', '--broad-cutoff', str(P))
    move_forward(WORK_DIR, macs_broad_p_suffix, ["*{}*".format(NAME), '*.bw', '*.bdg'], copy_only=True)
    process_macs2_logs(WORK_DIR + macs_broad_p_suffix)

Q = 0.01
# macs_q_suffix = "_macs_q{}".format(Q)
# if not os.path.exists(WORK_DIR + macs_q_suffix):
#     NAME = 'q{}'.format(Q)
#     run_bash("macs2.sh", WORK_DIR, GENOME, CHROM_SIZES, NAME, '-B', '-q', str(Q))
#     move_forward(WORK_DIR, macs_q_suffix, ["*{}*".format(NAME), '*.bw', '*.bdg'], copy_only=True)
#     process_macs2_logs(WORK_DIR + macs_q_suffix)

macs_broad_q_suffix = "_macs_broad_q{}".format(Q)
if not os.path.exists(WORK_DIR + macs_broad_q_suffix):
    NAME = 'broad_q{}'.format(Q)
    run_bash("macs2.sh", WORK_DIR, GENOME, CHROM_SIZES, NAME, '-B', '-q', str(Q), '--broad', '--broad-cutoff', str(Q))
    move_forward(WORK_DIR, macs_broad_q_suffix, ["*{}*".format(NAME), '*.bw', '*.bdg'], copy_only=True)
    process_macs2_logs(WORK_DIR + macs_broad_q_suffix)

# Custom fragment peak calling option
# FRAGMENT = 138
# macs_fragment_broad_suffix = "_macs_d{}_broad_q{}".format(FRAGMENT, Q)
# if not os.path.exists(WORK_DIR + macs_fragment_broad_suffix):
#     NAME = 'd{}_broad_q{}'.format(FRAGMENT, Q)
#     run_bash("macs2.sh", WORK_DIR, GENOME, CHROM_SIZES, NAME, '-B', '--broad', '--broad-cutoff', str(Q),
#              '--nomodel', '--shift', '0', '--extsize', str(FRAGMENT))
#     move_forward(WORK_DIR, macs_fragment_broad_suffix, ["*{}*".format(NAME), '*.bw', '*.bdg'], copy_only=True)
#     process_macs2_logs(WORK_DIR + macs_fragment_broad_suffix)

# # Batch macs14 with different peak calling procedures settings
# # P = 1e-5 is default for MACS14
# P = 0.00001
# macs14_suffix = '_macs14_p{}'.format(P)
# if not os.path.exists(WORK_DIR + macs14_suffix):
#     run_bash("macs14.sh", WORK_DIR, GENOME, str(P))
#     move_forward(WORK_DIR, macs14_suffix, ["*{}*".format(P)], copy_only=True)
#     process_macs2_logs(WORK_DIR + "_macs14_{}".format(P))

# Batch rseg
# rseg_suffix = '_rseg'
# if not os.path.exists(WORK_DIR + rseg_suffix):
#     run_bash("rseg.sh", WORK_DIR, GENOME, CHROM_SIZES)
#     move_forward(WORK_DIR, rseg_suffix, ["*-domains.bed", "*-scores.wig", "*-boundaries.bed",
#                                          "*-boundary-scores.wig", "*-counts.bed"], copy_only=True)
