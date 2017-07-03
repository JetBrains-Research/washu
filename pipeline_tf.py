#!/usr/bin/env python

# TODO:
# 1. GSM id map to SRX id
# 2. Download all (gsm, srx)
# +3. fastq-dump sras
# +4. fastqc + multiqc
# 5. bowties
# 6. macs2

import os
import itertools

import click
import pandas as pd

from reports.bowtie_logs import process_bowtie_logs
from pipeline_utils import *
from reports.peaks_logs import process_peaks_logs
from scripts.util import run_macs2


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
    GENOME = "hg19"
    INDEXES = os.path.join("/scratch/artyomov_lab_aging/Y20O20/chipseq/indexes",
                           GENOME)
    CHROM_SIZES = os.path.join(INDEXES, GENOME + ".chrom.sizes")

    # Data table
    data_table = pd.read_csv(data, sep="\t")

    # TODO: >>>>>>
    data_table = data_table.iloc[[0, 1], :]
    # TODO: <<<<<

    print("Data to process:")
    print(data_table)

    gsm2srxs = {}
    for r in data_table.itertuples():
        gsm2srxs[r.input] = r.input_srx.split(";")
        gsm2srxs[r.signal] = r.signal_srx.split(";")
    gsm_to_process = sorted(gsm2srxs.keys())
    print(gsm_to_process)

    # Make dirs:
    data_dirs = [os.path.join(out, gsmid) for gsmid in gsm_to_process]
    for data_dir in data_dirs:
        os.makedirs(data_dir, exist_ok=True)

    ##################
    # Pipeline start #
    ##################

    # Download SRA data:
    # 'rsync' here skips file if it already exist
    print("Downloading data...")
    srx_to_dir_list = []
    for gsmid in gsm_to_process:
        sra_dir = os.path.join(out, gsmid, "sra")
        srxs = gsm2srxs[gsmid]
        for srx in srxs:
            srx_to_dir_list.extend([srx, sra_dir])
    #TODO run_bash("geo_rsync.sh", *srx_to_dir_list)

    # Fastq-dump SRA data:
    #TODO run_bash("fastq_dump.sh", *data_dirs)

    # Prepare genome *.fa and Bowtie indexes
    print("Genomes and indices folder: ", INDEXES)
    run_bash("index_genome.sh", GENOME, INDEXES)
    run_bash("index_bowtie2.sh", GENOME, INDEXES)

    # Batch QC
    #TODO run_bash("fastqc.sh", *data_dirs)

    # Total multiqc:
    # Use -s options, otherwise tons of "SRRnnn" hard to distinguish
    # -s, --fullnames      Do not clean the sample names (leave as full
    #                      file name)
    # -f, --force          Overwrite any existing reports
    # -o, --outdir TEXT    Create report in the specified output directory.
    if len(data_dirs) > 1:
        run("multiqc", "-f", "-o", out, " ".join(data_dirs))

    # XXX: let's look in multiqc results, it shows that in several samples
    # it's better to trim first 5bp, so let's trim it in all samples for
    # simplicity:
    #
    #  * batch Bowtie with trim 5 first base pairs
    run_bash("index_bowtie2.sh", GENOME, INDEXES)
    run_bash("bowtie2.sh", GENOME, INDEXES, "5", *data_dirs)
    exit(-1)

    bams_dirs = []
    for data_dir in data_dirs:
        bams_dir = move_forward(data_dir, data_dir + "_bams",
                                ["*.bam", "*bowtie*.log"])
        bams_dirs.append(bams_dir)

        # multiqc is able to process Bowtie report
        run("multiqc", "-f", "-o", bams_dir, " ".join(data_dirs))

        # Create summary
        process_bowtie_logs(bams_dir)

    if len(data_dirs) > 1:
        run("multiqc", "-f", "-o", out, " ".join(data_dirs + bams_dirs))

    # Process insert size of BAM visualization
    # run_bash("fragments.sh", *bams_dirs)

    # Batch BigWig visualization
    # run_bash("bigwig.sh", WORK_DIR, CHROM_SIZES)
    # move_forward(WORK_DIR, WORK_DIR + "_bws", ["*.bw", "*.bdg", "*bw.log"],
    #              copy_only=True)
    #
    # # Batch RPKM visualization
    # run_bash("rpkm.sh", WORK_DIR)
    # move_forward(WORK_DIR, WORK_DIR + "_rpkms", ["*.bw", "*rpkm.log"],
    #              copy_only=True)
    #
    # # Remove duplicates
    # run_bash("remove_duplicates.sh", WORK_DIR)
    # move_forward(WORK_DIR, WORK_DIR + "_unique",
    #              ["*_unique*", "*_metrics.txt", "*duplicates.log"],
    #              copy_only=True)
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
    # folder = run_macs2(WORK_DIR, GENOME, CHROM_SIZES, 'broad_0.1', '--broad',
    #                    '--broad-cutoff', 0.1)
    # run_bash("../bed/macs2_filter_fdr.sh", folder,
    #          folder.replace('0.1', '0.05'), 0.1, 0.05, WORK_DIR)
    # run_bash("../bed/macs2_filter_fdr.sh", folder,
    #          folder.replace('0.1', '0.01'), 0.1, 0.01, WORK_DIR)
    #
    # # MACS2 Regular peak calling (https://github.com/taoliu/MACS) Q=0.01 in example
    # folder = run_macs2(WORK_DIR, GENOME, CHROM_SIZES, 'q0.1', '-q', 0.1)
    # run_bash("../bed/macs2_filter_fdr.sh", folder,
    #          folder.replace('0.1', '0.05'), 0.1, 0.05, WORK_DIR)
    # run_bash("../bed/macs2_filter_fdr.sh", folder,
    #          folder.replace('0.1', '0.01'), 0.1, 0.01, WORK_DIR)


if __name__ == '__main__':
    cli()
