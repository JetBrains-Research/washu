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
from logs.macs2_logs import process_macs2_logs
from logs.bowtie_logs import process_bowtie_logs
from pipeline_utils import *

parser = argparse.ArgumentParser(description='ULI ChIP-Seq data pipeline for WashU cluster')
parser.add_argument('path_to_directory', action=WritableDirectory, type=str,
                    help='Path to directory with data to run pipeline')
args = parser.parse_args()

# Configuration
WORK_DIR = args.path_to_directory
GENOME = "hg19"
INDEXES = os.path.join("/scratch/artyomov_lab_aging/indexes", GENOME)
CHROM_SIZES = os.path.join(INDEXES, GENOME + ".chrom.sizes")
READS = 15  # Subsampling to 15mln reads

print("Genomes and indices folder: ", INDEXES)
run_bash("index_genome.sh", GENOME, INDEXES)

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

# Batch BigWig visualization
run_bash("bigwig.sh", WORK_DIR, CHROM_SIZES)
move_forward(WORK_DIR, "_bws", ["*.bw", "*.bdg", "*bw.log"], copy_only=True)

# Batch subsampling
run_bash("subsample.sh", WORK_DIR, str(READS))
WORK_DIR = move_forward(WORK_DIR, "_{}mln".format(READS), ["*{}*".format(READS)])

# Batch BigWig visualization
run_bash("bigwig.sh", WORK_DIR, CHROM_SIZES)
move_forward(WORK_DIR, "_bws", ["*.bw", "*.bdg", "*bw.log"], copy_only=True)

# Batch macs with different peak calling procedures settings
for Q in [0.01, 0.001, 0.1]:
    run_bash("macs2.sh", WORK_DIR, GENOME, str(Q))
    move_forward(WORK_DIR, "_macs_{}".format(Q), ["*{}*".format(Q)], copy_only=True)
    process_macs2_logs(WORK_DIR + "_macs_{}".format(Q))

for Q in [0.01, 0.001, 0.1]:
    run_bash("macs2_broad.sh", WORK_DIR, GENOME, str(Q))
    move_forward(WORK_DIR, "_macs_broad_{}".format(Q), ["*{}*".format(Q)], copy_only=True)
    process_macs2_logs(WORK_DIR + "_macs_broad_{}".format(Q))
