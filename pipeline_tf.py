#!/usr/bin/env python

# TODO:
# 1. GSM id map to SRX id
# 2. Download all (gsm, srx)
# +3. fastq-dump sras
# +4. fastqc + multiqc
# 5. bowties
# 6. macs2

import os

import click
import pandas as pd

from pipeline_utils import *


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option('-d', '--data',
              required=True,
              type=click.Path(resolve_path=True, dir_okay=False, exists=True),
              help="Data *.tsv file")
@click.option('-o', '--out', default=".",
              type=click.Path(resolve_path=True, file_okay=False),
              help="Output dir (default: .)")
def cli(out, data):
    """ Download SRA data & call peaks

    \b
    VALUE is one of:
    """

    #################
    # Configuration #
    #################
    WORK_DIR = out
    GENOME = "hg19"
    INDEXES = os.path.join("/scratch/artyomov_lab_aging/Y20O20/chipseq/indexes", GENOME)
    CHROM_SIZES = os.path.join(INDEXES, GENOME + ".chrom.sizes")

    # Data table
    data_table = pd.read_csv(data, sep="\t")
    data_table = data_table.iloc[[0, 1], :]

    print("Data to process:")
    print(data_table)

    gsm2srx = {}
    for r in data_table.itertuples():
        gsm2srx[r.input] = r.input_srx
        gsm2srx[r.signal] = r.signal_srx
    gsm_to_process = sorted(gsm2srx.keys())
    print(gsm_to_process)

    # Make dirs:
    data_dirs = [os.path.join(out, gsmid) for gsmid in gsm_to_process]
    for d in data_dirs:
        os.makedirs(d, exist_ok=True)

    ##################
    # Pipeline start #
    ##################

    # Download SRA data:
    # 'rsync' here skips file if it already exist
    for gsmid in gsm_to_process:

        sra_dir = os.path.join(out, gsmid, "sra")
        os.makedirs(sra_dir, exist_ok=True)

        print("Downloading data to {} ...".format(sra_dir))
        srxid = gsm2srx[gsmid]
        url = "rsync://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/" \
              "sra/SRX/{}/{}".format(srxid[0:6], srxid)

        # Options:
        #   -a, --archive            archive mode; equals -rlptgoD (no -H,-A,-X)
        #   -z, --compress           compress file data during the transfer
        #   -v, --verbose            increase verbosity
        #   -i,  --itemize-changes   output a change-summary for all updates
        #   -t, --times              preserve modification times
        #   —partial (or -P = —progress —partial): enable partial transmission,
        #                            but not for local fs
        #   -r, --recursive          recurse into directories
        #   -h, --human-readable     output numbers in a human-readable format
        run("rsync -azvvit --partial-dir=.rsync-partial --human-readable"
            " --progress {} {}".format(url, sra_dir))

    # Fastq-dump SRA data:
    run_bash("fastq_dump.sh", *data_dirs)

    # Prepare genome *.fa and Bowtie indexes
    print("Genomes and indices folder: ", INDEXES)
    run_bash("index_genome.sh", GENOME, INDEXES)
    run_bash("index_bowtie.sh", GENOME, INDEXES)

    # Batch QC
    run_bash("fastqc.sh", *data_dirs)

    # Total multiqc:
    if len(data_dirs) > 1:
        run("multiqc", "-f", "-o", WORK_DIR, " ".join(data_dirs))

    # # Batch Bowtie with trim 5 first base pairs
    # run_bash("index_bowtie.sh", GENOME, INDEXES)
    # run_bash("bowtie.sh", WORK_DIR, GENOME, INDEXES, "5")
    # WORK_DIR = move_forward(WORK_DIR, WORK_DIR + "_bams", ["*.bam", "*bowtie*.log"])
    # # multiqc is able to process Bowtie report
    # subprocess.run("multiqc " + WORK_DIR, shell=True)
    # # Create summary
    # process_bowtie_logs(WORK_DIR)
    #
    # # Process insert size of BAM visualization
    # run_bash("fragments.sh", WORK_DIR)
    #
    # # Batch BigWig visualization
    # run_bash("bigwig.sh", WORK_DIR, CHROM_SIZES)
    # move_forward(WORK_DIR, WORK_DIR + "_bws", ["*.bw", "*.bdg", "*bw.log"], copy_only=True)
    #
    # # Batch RPKM visualization
    # run_bash("rpkm.sh", WORK_DIR)
    # move_forward(WORK_DIR, WORK_DIR + "_rpkms", ["*.bw", "*rpkm.log"], copy_only=True)
    #
    # # Remove duplicates
    # run_bash("remove_duplicates.sh", WORK_DIR)
    # move_forward(WORK_DIR, WORK_DIR + "_unique", ["*_unique*", "*_metrics.txt", "*duplicates.log"], copy_only=True)
    #
    # # Batch subsampling to 15mln reads
    # # READS = 15
    # # run_bash("subsample.sh", WORK_DIR, str(READS))
    # # WORK_DIR = move_forward(WORK_DIR, WORK_DIR + "_{}mln".format(READS), ["*{}*".format(READS)])
    # #
    # # # Batch BigWig visualization
    # # run_bash("bigwig.sh", WORK_DIR, CHROM_SIZES)
    # # move_forward(WORK_DIR, WORK_DIR + "_bws", ["*.bw", "*.bdg", "*bw.log"], copy_only=True)
    #
    # ########################
    # # Peak calling section #
    # ########################
    #
    # # MACS2 Broad peak calling (https://github.com/taoliu/MACS) Q=0.1 in example
    # folder = run_macs2(WORK_DIR, GENOME, CHROM_SIZES, 'broad_0.1', '--broad', '--broad-cutoff', 0.1)
    # run_bash("../bed/macs2_filter_fdr.sh", folder, folder.replace('0.1', '0.05'), 0.1, 0.05, WORK_DIR)
    # run_bash("../bed/macs2_filter_fdr.sh", folder, folder.replace('0.1', '0.01'), 0.1, 0.01, WORK_DIR)
    #
    # # MACS2 Regular peak calling (https://github.com/taoliu/MACS) Q=0.01 in example
    # folder = run_macs2(WORK_DIR, GENOME, CHROM_SIZES, 'q0.1', '-q', 0.1)
    # run_bash("../bed/macs2_filter_fdr.sh", folder, folder.replace('0.1', '0.05'), 0.1, 0.05, WORK_DIR)
    # run_bash("../bed/macs2_filter_fdr.sh", folder, folder.replace('0.1', '0.01'), 0.1, 0.01, WORK_DIR)

if __name__ == '__main__':
    cli()
