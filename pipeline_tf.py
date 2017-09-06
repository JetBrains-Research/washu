#!/usr/bin/env python

# TODO: Encode pipeline suggest using bowtie1 + SPP aligner for TF data

import sys

import click
import pandas as pd

from pipeline_utils import *
from reports.bowtie_logs import process_bowtie_logs
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
    """ For given descriptor table download SRA data & call peaks

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

    # TODO: clean-up when pipeline will be finished
    # TODO >>>>>>
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

    run_bash("geo_rsync.sh", *srx_to_dir_list)

    # Fastq-dump SRA data:
    run_bash("fastq_dump.sh", *data_dirs)

    # Prepare genome *.fa and Bowtie indexes
    print("Genomes and indices folder: ", INDEXES)
    run_bash("index_genome.sh", GENOME, INDEXES)
    run_bash("index_bowtie2.sh", GENOME, INDEXES)
    #run_bash("index_bowtie.sh", GENOME, INDEXES)

    # Batch QC
    run_bash("fastqc.sh", *data_dirs)

    # Total multiqc:
    # Use -s options, otherwise tons of "SRRnnn" hard to distinguish
    # -s, --fullnames      Do not clean the sample names (leave as full
    #                      file name)
    # -f, --force          Overwrite any existing reports
    # -o, --outdir TEXT    Create report in the specified output directory.
    # if len(data_dirs) > 1:
    #     run("multiqc", "-f", "-o", out, " ".join(data_dirs))
    #
    # XXX: let's look in multiqc results, it shows that in several samples
    # it's better to trim first 5bp, so let's trim it in all samples for
    # simplicity:
    #
    #  * batch Bowtie with trim 5 first base pairs
    run_bash("bowtie2.sh", GENOME, INDEXES, "5", *data_dirs)

    # Merge TF SRR*.bam files to one
    run_bash("samtools_merge.sh", GENOME, *data_dirs)

    bams_dirs = []
    for data_dir in data_dirs:
        bams_dir = data_dir + "_bams"
        bams_dirs.append(bams_dir)
        move_forward(data_dir, bams_dir, ["*.bam", "*bowtie*.log"])

        # multiqc is able to process Bowtie report
        run("multiqc", "-f", "-o", bams_dir, " ".join(bams_dirs))

        # Create summary
        process_bowtie_logs(bams_dir)

    if len(data_dirs) > 1:
        run("multiqc", "-f", "-o", out, " ".join(data_dirs + bams_dirs))

    # XXX: doesn't work for some reason, "filter by -f66" returns nothing
    # Process insert size of BAM visualization
    # run_bash("fragments.sh", *bams_dirs)

    # Batch BigWig visualization
    run_bash("bigwig.sh", CHROM_SIZES, *bams_dirs)
    for bams_dir in bams_dirs:
        move_forward(bams_dir, bams_dir + "_bws", ["*.bw", "*.bdg", "*bw.log"],
                     copy_only=True)

    # Batch RPKM visualization
    run_bash("rpkm.sh", *bams_dirs)
    for bams_dir in bams_dirs:
        move_forward(bams_dir, bams_dir + "_rpkms", ["*.bw", "*rpkm.log"],
                     copy_only=True)

    # Remove duplicates
    run_bash("remove_duplicates.sh", *bams_dirs)
    for bams_dir in bams_dirs:
        move_forward(bams_dir, bams_dir + "_unique",
                     ["*_unique*", "*_metrics.txt", "*duplicates.log"],
                     copy_only=True)

    # Call PEAKS:
    files_to_cleanup = []
    try:
        for r in data_table.itertuples():
            gsmid_signal = r.signal
            gsmid_input = r.input

            bams_dir_signal = os.path.join(out, gsmid_signal + "_bams")
            bams_dir_input = os.path.join(out, gsmid_input + "_bams")

            # Find all input *.bam and *.bam.bai
            input_files = [f for f in os.listdir(bams_dir_input)
                           if f.endswith(".bam") or f.endswith(".bam.bai")]
            # Create symlink to link
            for f in input_files:
                f_link = os.path.join(bams_dir_signal,
                                      f.replace(".bam", "_input.bam"))
                run("ln", "-s", os.path.join(bams_dir_input, f), f_link)
                files_to_cleanup.append(f_link)

            ########################
            # Peak calling section #
            ########################
            # TODO: Vectorize QSUB MACS2 script:

            # MACS2 Broad peak calling (https://github.com/taoliu/MACS) Q=0.1
            #  in example
            peaks_dirs = run_macs2(
                GENOME, CHROM_SIZES,
                'broad_0.1', '--broad', '--broad-cutoff', 0.1,
                work_dirs=[bams_dir_signal])

            # Bedtools is necessary for filter script
            subprocess.run('module load bedtools2', shell=True)
            for peaks_dir in peaks_dirs:
                run_bash('../bed/macs2_filter_fdr.sh', peaks_dir,
                         peaks_dir.replace('0.1', '0.05'), 0.1, 0.05,
                         bams_dir_signal)
                run_bash('../bed/macs2_filter_fdr.sh', peaks_dir,
                         peaks_dir.replace('0.1', '0.01'), 0.1, 0.01,
                         bams_dir_signal)

            # # MACS2 Regular peak calling (https://github.com/taoliu/MACS)
            # # Q=0.01 in example
            peaks_dirs = run_macs2(
                GENOME, CHROM_SIZES, 'q0.1', '-q', 0.1,
                work_dirs=[bams_dir_signal])
            for peaks_dir in peaks_dirs:
                run_bash("../bed/macs2_filter_fdr.sh", peaks_dir,
                         peaks_dir.replace('0.1', '0.05'), 0.1, 0.05,
                         bams_dir_signal)
                run_bash("../bed/macs2_filter_fdr.sh", peaks_dir,
                         peaks_dir.replace('0.1', '0.01'), 0.1, 0.01,
                         bams_dir_signal)
    finally:
        for f in files_to_cleanup:
            print("Cleanup:")
            try:
                os.remove(f)
                print("* deleted: ", f)
            except OSError:
                print("Error while deleting '{}'".format(f), sys.exc_info()[0])

if __name__ == '__main__':
    cli()
