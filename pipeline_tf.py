#!/usr/bin/env python

# TODO: Encode pipeline suggest using bowtie1 + SPP aligner for TF data
# TODO: table 'folder' value not used

import sys

import pandas as pd

from pipeline_utils import *
from reports.bowtie_logs import process_bowtie_logs
from scripts.util import run_macs2


def cli():
    parser = argparse.ArgumentParser(
        description="For given descriptor table download SRA data & call peaks"
    )
    parser.add_argument("data",
                        # required=True,
                        type=file_path_type(ext="tsv"),
                        help="Data *.tsv file")
    parser.add_argument('-o', '--out', default='.',
                        type=file_path_type(True, False),
                        help="Output dir (default: .)")

    args = parser.parse_args()

    run_pipeline(args.out, args.data)


def run_pipeline(out, data):
    #################
    # Configuration #
    #################
    genome = "hg19"
    indexes \
        = os.path.join("/scratch/artyomov_lab_aging/Y20O20/chipseq/indexes",
                       genome)
    chrom_sizes = os.path.join(indexes, genome + ".chrom.sizes")
    picard_tools = os.path.join("~", "picard.jar")

    # Data table
    data_table = pd.read_csv(data, sep="\t")

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

    run_bash("scripts/geo_rsync.sh", *srx_to_dir_list)

    # Fastq-dump SRA data:
    run_bash("parallel/fastq_dump.sh", *data_dirs)

    # Prepare genome *.fa and Bowtie indexes
    print("Genomes and indices folder: ", indexes)
    run_bash("parallel/index_genome.sh", genome, indexes)
    run_bash("parallel/index_bowtie2.sh", genome, indexes)
    # run_bash("index_bowtie.sh", genome, indexes)

    # Batch QC
    run_bash("parallel/fastqc.sh", *data_dirs)

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

    # Alignment step:
    def process_sra(sra_dirs):
        #  * batch Bowtie with trim 5 first base pairs
        run_bash("parallel/bowtie2.sh", genome, indexes, "5", *sra_dirs)

        # Merge TF SRR*.bam files to one
        run_bash("parallel/samtools_merge.sh", genome, *sra_dirs)

    bams_dirs = process_dirs(data_dirs, "_bams", ["*.bam", "*bowtie*.log"],
                             process_sra)
    for bams_dir in bams_dirs:
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
    process_dirs(bams_dirs, "_bws", ["*.bw", "*.bdg", "*bw.log"],
                 lambda dirs: run_bash("parallel/bigwig.sh", chrom_sizes,
                                       *dirs))

    # Batch RPKM visualization
    process_dirs(bams_dirs, "_rpkms", ["*.bw", "*rpkm.log"],
                 lambda dirs: run_bash("parallel/rpkm.sh", *dirs))

    # Remove duplicates
    process_dirs(bams_dirs, "_unique",
                 ["*_unique*", "*_metrics.txt", "*duplicates.log"],
                 lambda dirs: run_bash("parallel/remove_duplicates.sh",
                                       picard_tools, *dirs))

    # Call PEAKS:
    files_to_cleanup = []
    try:
        # let's link signal bams with corresponding input:
        bams_dirs_for_peakcalling = []
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

            bams_dirs_for_peakcalling.append(bams_dir_signal)

        ########################
        # Peak calling section #
        ########################
        # Bedtools is necessary for filter script
        subprocess.run('module load bedtools2', shell=True)

        # MACS2 Broad peak calling (https://github.com/taoliu/MACS) Q=0.1
        #  in example
        peaks_dirs = run_macs2(genome, chrom_sizes,
                               'broad_0.1', '--broad', '--broad-cutoff', 0.1,
                               work_dirs=bams_dirs_for_peakcalling)
        for bams_dir_signal, peaks_dir in zip(bams_dirs_for_peakcalling,
                                              peaks_dirs):
            run_bash('bed/macs2_filter_fdr.sh', peaks_dir,
                     peaks_dir.replace('0.1', '0.05'), 0.1, 0.05,
                     bams_dir_signal)
            run_bash('bed/macs2_filter_fdr.sh', peaks_dir,
                     peaks_dir.replace('0.1', '0.01'), 0.1, 0.01,
                     bams_dir_signal)

        # # MACS2 Regular peak calling (https://github.com/taoliu/MACS)
        # # Q=0.01 in example
        peaks_dirs = run_macs2(genome, chrom_sizes, 'q0.1', '-q', 0.1,
                               work_dirs=bams_dirs_for_peakcalling)
        for bams_dir_signal, peaks_dir in zip(bams_dirs_for_peakcalling,
                                              peaks_dirs):
            run_bash("bed/macs2_filter_fdr.sh", peaks_dir,
                     peaks_dir.replace('0.1', '0.05'), 0.1, 0.05,
                     bams_dir_signal)
            run_bash("bed/macs2_filter_fdr.sh", peaks_dir,
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


def process_dirs(dirs, suffix, what_to_move, processor_fun):
    # filter already processed dirs:
    dirs_to_process = []
    result_dirs = []
    for data_dir in dirs:
        res_dir = data_dir + suffix
        result_dirs.append(res_dir)
        if os.path.exists(res_dir):
            print("[Skipped] Directory already exists:", res_dir)
            continue
        dirs_to_process.append((data_dir, res_dir))

    if dirs_to_process:
        # process dirs:
        processor_fun(list(zip(*dirs_to_process))[0])

        # move results:
        for data_dir, res_dir in dirs_to_process:
            move_forward(data_dir, res_dir, what_to_move, copy_only=True)

    return result_dirs


if __name__ == '__main__':
    cli()
