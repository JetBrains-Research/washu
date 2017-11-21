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
from reports.bowtie_logs import process_bowtie_logs
from reports.peaks_logs import process_peaks_logs
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
WORK_DIR = args.path_to_directory
GENOME = args.genome
INDEXES = os.path.join(args.path_to_indexes, GENOME)

print("WORK_DIR:", WORK_DIR)
print("GENOME:", GENOME)
print("INDEXES:", INDEXES)

##################
# Pipeline start #
##################
CHROM_SIZES = os.path.join(INDEXES, GENOME + ".chrom.sizes")
PICARD_TOOLS = os.path.join("~", "picard.jar")
run_bash("parallel/index_genome.sh", GENOME, INDEXES)

# Batch QC
run_bash("parallel/fastqc.sh", WORK_DIR)

# Batch Bowtie with trim 5 first base pairs
run_bash("parallel/index_bowtie.sh", GENOME, INDEXES)
run_bash("parallel/bowtie.sh", GENOME, INDEXES, "5", WORK_DIR)
move_forward(WORK_DIR, WORK_DIR + "_bams", ["*.bam", "*bowtie*.log"])
WORK_DIR = WORK_DIR + "_bams"
os.chdir(WORK_DIR)

# multiqc is able to process Bowtie report
subprocess.run("multiqc " + WORK_DIR, shell=True)
# Create summary
process_bowtie_logs(WORK_DIR)

# Batch BigWig visualization
run_bash("parallel/bigwig.sh", CHROM_SIZES, WORK_DIR)
move_forward(WORK_DIR, WORK_DIR + "_bws", ["*.bw", "*.bdg", "*bw.log"])

# Remove duplicates
run_bash("parallel/remove_duplicates.sh", PICARD_TOOLS, WORK_DIR)
move_forward(WORK_DIR, WORK_DIR + "_unique",
             ["*_unique*", "*_metrics.txt", "*duplicates.log"])

# Tags BW visualization
run_bash("parallel/tags_bigwig.sh", CHROM_SIZES, 150, WORK_DIR + "_unique")
move_forward(WORK_DIR + "_unique", WORK_DIR + "_unique_tags_bws", ["*bw*"])

# Batch RPKM visualization
run_bash("parallel/rpkm.sh", WORK_DIR)
move_forward(WORK_DIR, WORK_DIR + "_rpkms", ["*.bw", "*rpkm.log"])

# Batch subsampling to 15mln reads
# READS = 15
# run_bash("subsample.sh", WORK_DIR, str(READS))
# WORK_DIR = move_forward(WORK_DIR, WORK_DIR + "_{}mln".format(READS),
#                         ["*{}*".format(READS)])
#
# # Batch BigWig visualization
# run_bash("bigwig.sh", CHROM_SIZES, WORK_DIR)
# move_forward(WORK_DIR, WORK_DIR + "_bws", ["*.bw", "*.bdg", "*bw.log"],
#              chdir=False)

########################
# Peak calling section #
########################

# MACS2 Broad peak calling (https://github.com/taoliu/MACS) Q=0.1 in example
folder = run_macs2(GENOME, CHROM_SIZES,
                   'broad_0.1', '--broad', '--broad-cutoff', 0.1,
                   work_dirs=[WORK_DIR])[0]
peaks_folder = folder.replace('0.1', '0.05')
run_bash('bed/macs2_filter_fdr.sh', folder, peaks_folder, 0.1, 0.05, WORK_DIR)
run_bash('parallel/peaks_frip.sh', peaks_folder, WORK_DIR)

peaks_folder = folder.replace('0.1', '0.01')
run_bash('bed/macs2_filter_fdr.sh', folder, peaks_folder, 0.1, 0.01, WORK_DIR)
run_bash('parallel/peaks_frip.sh', peaks_folder, WORK_DIR)

# MACS2 Regular peak calling (https://github.com/taoliu/MACS) Q=0.01 in example
folder = run_macs2(GENOME, CHROM_SIZES, 'q0.1', '-q', 0.1,
                   work_dirs=[WORK_DIR])[0]
peaks_folder = folder.replace('0.1', '0.05')
run_bash('bed/macs2_filter_fdr.sh', folder, peaks_folder, 0.1, 0.05, WORK_DIR)
run_bash('parallel/peaks_frip.sh', peaks_folder, WORK_DIR)

peaks_folder = folder.replace('0.1', '0.01')
run_bash('bed/macs2_filter_fdr.sh', folder, peaks_folder, 0.1, 0.01, WORK_DIR)
run_bash('parallel/peaks_frip.sh', peaks_folder, WORK_DIR)

# MACS1.4 P=1e-5 is default
# P = 0.00001
# NAME = '14_p{}'.format(P)
# FOLDER = '{}_macs_{}'.format(WORK_DIR, NAME)
# print(FOLDER)
# if not os.path.exists(FOLDER):
#     run_bash("macs14.sh", WORK_DIR, GENOME, str(P))
#     move_forward(WORK_DIR, FOLDER, ['*{}*'.format(NAME), '*rip.csv'],
#                  chdir=False)
#     process_macs2_logs(FOLDER)

# Batch RSEG
rseg_suffix = '_rseg'
if not os.path.exists(WORK_DIR + rseg_suffix):
    run_bash("parallel/rseg.sh", WORK_DIR, GENOME, CHROM_SIZES)
    move_forward(WORK_DIR, WORK_DIR + rseg_suffix,
                 ['*domains*', '*rseg*', '*.bam.bed', 'deadzones*',
                  '*_chrom_sizes.bed', '*rip.csv'])
    process_peaks_logs(WORK_DIR + rseg_suffix)

# Batch SICER
sicer_suffix = '_sicer'
if not os.path.exists(WORK_DIR + sicer_suffix):
    # <work_dir> <genome> <chrom.sizes> <FDR> [window size (bp)] [fragment size] [gap size (bp)] [batch] # nopep8
    run_bash("parallel/sicer.sh", WORK_DIR, GENOME, CHROM_SIZES, "0.01", "200", "150", "0")
    move_forward(WORK_DIR, WORK_DIR + sicer_suffix,
                 ['*sicer.log', '*island.bed', '*rip.csv'])
    process_peaks_logs(WORK_DIR + sicer_suffix)
