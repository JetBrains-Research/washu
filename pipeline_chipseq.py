#!/usr/bin/env python

"""
This is a chip-seq technical pipeline.

Usage:
  * FastQC
  * Visualization
  * Decide whether trimming or subsampling is required
  * Modify inplace copy of this pipeline
  * Launch pipeline

Conventions:
This pipeline uses folder naming as a steps, i.e. next step appends _suffix for the working folder,
stores results in new folder and change working folder if necessary.

NOTE: python3 required
> source activate py3.5

author oleg.shpynov@jetbrains.com
"""
from pipeline_utils import *
from scripts.util import run_macs2

parser = argparse.ArgumentParser(description='ULI ChIP-Seq data pipeline')
parser.add_argument('path_to_directory', action=WritableDirectory, type=str,
                    help='Path to directory with data to run pipeline')
parser.add_argument('path_to_indexes', action=WritableDirectory, type=str,
                    help='Path to indexes')
parser.add_argument('genome', type=str,
                    help='Genome')
args = parser.parse_args()

#################
# Configuration #
#################
WORK_DIR = args.path_to_directory.rstrip('/')
GENOME = args.genome
INDEXES = os.path.join(args.path_to_indexes, GENOME)
CHROM_SIZES = os.path.join(INDEXES, GENOME + ".chrom.sizes")

print("WORK_DIR:", WORK_DIR)
print("GENOME:", GENOME)
print("INDEXES:", INDEXES)
print("CHROM_SIZES:", CHROM_SIZES)

PICARD_TOOLS = os.path.expanduser("~/picard.jar")
print("PICARD_TOOLS:", PICARD_TOOLS)
PHANTOMPEAKQUALTOOLS = os.path.expanduser("~/phantompeakqualtools")
print("PHANTOMPEAKQUALTOOLS:", PHANTOMPEAKQUALTOOLS)
SPAN = os.path.expanduser("~/span.jar")
print("SPAN:", SPAN)

##################
# Pipeline start #
##################
run_bash("parallel/index_genome.sh", GENOME, INDEXES)

# Batch QC
run_bash("parallel/fastqc.sh", WORK_DIR)
# multiqc
subprocess.run("multiqc -f " + WORK_DIR, shell=True)

# Batch Bowtie with trim 5 first base pairs
run_bash("parallel/index_bowtie.sh", GENOME, INDEXES)
run_bash("parallel/bowtie.sh", GENOME, INDEXES, "5", WORK_DIR)
move_forward(WORK_DIR, WORK_DIR + "_bams", ["*.bam", "*bowtie*.log"])
# multiqc
subprocess.run("multiqc -f " + WORK_DIR + "_bams", shell=True)

WORK_DIR = WORK_DIR + "_bams"
os.chdir(WORK_DIR)

# Batch BigWig visualization
run_bash("parallel/bigwig.sh", CHROM_SIZES, WORK_DIR)
move_forward(WORK_DIR, WORK_DIR + "_bws", ["*.bw", "*.bdg", "*bw.log"])

# QC PBC/NRF + PhantomPeakQualTools metrics for BAMs
run_bash("parallel/bam_qc.sh", PHANTOMPEAKQUALTOOLS, WORK_DIR)
move_forward(WORK_DIR, WORK_DIR + "/qc", ["*.pdf", "*phantom.tsv", "*pbc_nrf.tsv", "*bam_qc*"])

# Remove duplicates
run_bash("parallel/remove_duplicates.sh", PICARD_TOOLS, WORK_DIR)
move_forward(WORK_DIR, WORK_DIR + "_unique",
             ["*_unique*", "*_metrics.txt", "*duplicates.log"])

# Batch subsampling to 15mln reads
# READS = 15
# run_bash("subsample.sh", WORK_DIR, str(READS))
# WORK_DIR = move_forward(WORK_DIR, WORK_DIR + "_{}mln".format(READS),
#                         ["*{}*".format(READS)])

########################
# Peak calling section #
########################

# MACS2 Broad peak calling. Default: Q=0.1
run_macs2(GENOME, CHROM_SIZES, 'broad_0.1', '--broad', '--broad-cutoff',
          0.1, work_dirs=[WORK_DIR])

# MACS2 Regular peak calling. Default: Q=0.05
run_macs2(GENOME, CHROM_SIZES, 'q0.05', '-q', 0.05, work_dirs=[WORK_DIR])

# MACS1.4 P=1e-5 is default
# P = 0.00001
# NAME = '14_p{}'.format(P)
# FOLDER = '{}_macs_{}'.format(WORK_DIR, NAME)
# print(FOLDER)
# if not os.path.exists(FOLDER):
#     run_bash("macs14.sh", WORK_DIR, GENOME, str(P))
#     move_forward(WORK_DIR, FOLDER, ['*{}*'.format(NAME), '*rip.csv'],
#                  chdir=False)

# Batch RSEG
rseg_suffix = '_rseg'
if not os.path.exists(WORK_DIR + rseg_suffix):
    run_bash("parallel/rseg.sh", WORK_DIR, GENOME, CHROM_SIZES)
    move_forward(WORK_DIR, WORK_DIR + rseg_suffix,
                 ['*domains*', '*rseg*', '*.bam.bed', 'deadzones*',
                  '*_chrom_sizes.bed', '*rip.csv'])
    # multiqc
    subprocess.run("multiqc " + WORK_DIR + rseg_suffix, shell=True)

# Batch SICER
sicer_suffix = '_sicer'
if not os.path.exists(WORK_DIR + sicer_suffix):
    # <work_dir> <genome> <chrom.sizes> <FDR> [window size (bp)] [fragment size] [gap size (bp)] # nopep8
    run_bash("parallel/sicer.sh", WORK_DIR, GENOME, CHROM_SIZES, "0.01", "200", "150", "0")
    move_forward(WORK_DIR, WORK_DIR + sicer_suffix,
                 ['*sicer.log', '*.bed', '*rip.csv'])
    # multiqc
    subprocess.run("multiqc " + WORK_DIR + sicer_suffix, shell=True)

# Batch SPAN
span_suffix = '_span'
if not os.path.exists(WORK_DIR + span_suffix):
    os.mkdir(WORK_DIR + span_suffix)
    #  <SPAN_JAR_PATH> <WORK_DIR> <GENOME> <CHROM_SIZES> <Q> [<OUTPUT_DIR> [<GAP>]]
    run_bash("parallel/span.sh", SPAN, WORK_DIR, GENOME, CHROM_SIZES, "0.01",
             WORK_DIR + span_suffix, 5)
    move_forward(WORK_DIR, WORK_DIR + span_suffix, ['*_peaks.bed', '*span*.log'])
