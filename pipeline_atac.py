import os
import os.path
import subprocess
import shutil

from bowtie_logs import process as bowtie_process
from macs2_logs import process as macs2_process
import argparse


class WritableDirectory(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if not os.path.isdir(values):
            raise argparse.ArgumentTypeError("{0} is not a valid path".format(values))
        if os.access(values, os.W_OK):
            setattr(namespace, self.dest ,values)
        else:
            raise argparse.ArgumentTypeError("{0} is not a writable directory".format(values))

parser = argparse.ArgumentParser(description='ATAC-seq data pipeline for WashU cluster')
parser.add_argument('path_to_directory', action=WritableDirectory, type=str,
                    help='Path to directory with data to run pipeline')

args = parser.parse_args()

SCRIPTS_PATH = "~/work/washu/scripts/"


def run_bash(script_file, *args):
    command = " ".join(["bash",
                        os.path.join(SCRIPTS_PATH, script_file),
                        *args])
    print(command)
    subprocess.run(command, shell=True)


def move_results(folder, what_to_move):
    os.mkdir(folder)
    for pattern in what_to_move:
        shutil.move(pattern, folder)


# Configuration
WORK_DIR = args.path_to_directory
GENOME = "hg19" # Nothing to compare with aligned on hg38
INDEXES = WORK_DIR + "/../" + GENOME
CHROMSIZES = os.path.join(INDEXES, GENOME + ".chrom.sizes")


print("Genomes and indices folder: " + INDEXES)
run_bash("genome_indices.sh", GENOME, INDEXES)
os.chdir(WORK_DIR)

# Batch QC
run_bash("fastqc.sh", WORK_DIR)

run_bash("bowtie2.sh", WORK_DIR, GENOME, INDEXES)
BAMS = WORK_DIR + "_bams"
move_results(BAMS, ["*.bam", "*bowtie*.log"])
os.chdir(BAMS)

# Batch fragments
run_bash("fragments.sh", WORK_DIR)
move_results(WORK_DIR + "_fragments", ["InsertSizeMetrics.txt", "fragments.png"])

# Create summary
bowtie_process(BAMS)
WORK_DIR = os.getcwd()
print("Working directory: " + WORK_DIR)

# MultiQC is able to process Bowtie report
print("Processing multiqc")
subprocess.run("multiqc " + WORK_DIR, shell=True)

# Batch BigWig visualization
run_bash("bigwig.sh", WORK_DIR, CHROMSIZES)
BWS = WORK_DIR + "_bws"
move_results(BWS, ["*.bw", "*.bdg", "*bw.log"])

# Batch subsampling
READS = 15
print("Subsampling to {}mln".format(READS))
run_bash("sabsample.sh", WORK_DIR, READS)

# Move results and CD
SUBSAMPLED = "{}_{}mln".format(WORK_DIR, str(READS))
move_results(SUBSAMPLED, "*{}*".format(str(READS)))
os.chdir(SUBSAMPLED)
WORK_DIR = os.getcwd()
print("Working directory: " + WORK_DIR)

# Batch BigWig visualization
run_bash("bigwig.sh", WORK_DIR, CHROMSIZES)
BWS = WORK_DIR + "_bws"
move_results(BWS, ["*.bw", "*.bdg", "*bw.log"])

# Batch macs with different peak calling procedures settings
QS = (0.001, 0.01, 0.1)
for Q in QS:
    run_bash("macs2.sh", WORK_DIR, GENOME, Q, CHROMSIZES)
    peaks = "{}_macs_{}".format(WORK_DIR, str(Q))
    move_results(peaks, "*{}*".format(str(Q)))
    macs2_process(peaks)

# QS = (0.001, 0.01, 0.1)
# for Q in QS:
#     run_bash("macs2_broad.sh", WORK_DIR, GENOME, Q, CHROMSIZES)
#     peaks = "{}_macs_broad_{}".format(WORK_DIR, str(Q))
#     move_results(peaks, "*{}*".format(str(Q)))
#     macs2_process(peaks)
#
# print("Done")