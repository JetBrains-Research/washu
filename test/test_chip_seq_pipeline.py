import unittest
from subprocess import call
import os
import shutil
import glob


def run_pipeline():
    os.chdir(os.path.expanduser("~"))

    call(["tar", "xfz", "chip-seq-pipeline-test-data.tar.gz"])

    call(["mv", "chip-seq-pipeline-test-data", "test_data"])

    shutil.copytree("./test_data/fastq/", "./H3K4me3")

    os.makedirs("/scratch/artyomov_lab_aging/Y20O20/chipseq/indexes")

    shutil.copytree("./test_data/index/hg19", "/scratch/artyomov_lab_aging/Y20O20/chipseq/indexes/hg19")

    os.mkdir("./H3K4me3_bams")

    shutil.copy("./test_data/deadzones-k36-hg19.bed", "./H3K4me3_bams/deadzones-k36-hg19.bed")

    call(["python", "/washu/pipeline_chipseq.py", "/root/H3K4me3"])


class ChipSeqPipeline(unittest.TestCase):
    def check_files_not_empty(self, pattern, expected_files_number=None):
        files = glob.glob(pattern)
        self.assertEqual(expected_files_number, len(files),
                         "Expected {} files for pattern '{}', but {} found.".format(expected_files_number, pattern,
                                                                                    len(files)))

        for f in files:
            self.assertNotEqual(0, os.path.getsize(f), "File {} is empty.".format(f))

    def check_files(self):
        os.chdir(os.path.expanduser("~"))
        self.check_files_not_empty("./H3K4me3_bams/*.bam", 5)
        self.check_files_not_empty("./H3K4me3_bams/*.bai", 5)

        self.check_files_not_empty("./H3K4me3_bams_macs2_broad_0.01/*.broadPeak", 4)
        self.check_files_not_empty("./H3K4me3_bams_macs2_broad_0.01/*.broadPeak_rip.csv", 4)

        self.check_files_not_empty("./H3K4me3_bams_macs2_broad_0.05/*.broadPeak", 4)
        self.check_files_not_empty("./H3K4me3_bams_macs2_broad_0.05/*.broadPeak_rip.csv", 4)

        self.check_files_not_empty("./H3K4me3_bams_macs2_broad_0.1/*.broadPeak", 4)
        self.check_files_not_empty("./H3K4me3_bams_macs2_broad_0.1/*.broadPeak_rip.csv", 4)

        self.check_files_not_empty("./H3K4me3_bams_macs2_q0.01/*.narrowPeak", 4)
        self.check_files_not_empty("./H3K4me3_bams_macs2_q0.01/*.narrowPeak_rip.csv", 4)

        self.check_files_not_empty("./H3K4me3_bams_macs2_q0.05/*.narrowPeak", 4)
        self.check_files_not_empty("./H3K4me3_bams_macs2_q0.05/*.narrowPeak_rip.csv", 4)

        self.check_files_not_empty("./H3K4me3_bams_macs2_q0.1/*.narrowPeak", 4)
        self.check_files_not_empty("./H3K4me3_bams_macs2_q0.1/*.narrowPeak_rip.csv", 4)

        self.check_files_not_empty("./H3K4me3_bams_rpkms/*.bw", 5)

        self.check_files_not_empty("./H3K4me3_bams_rseg/*_domains.bed", 4)

        self.check_files_not_empty("./H3K4me3_bams_sicer_0.01/*-island.bed", 4)
        self.check_files_not_empty("./H3K4me3_bams_sicer_0.01/*-island.bed_rip.csv", 4)

        self.check_files_not_empty("./H3K4me3_bams_unique/*.bam", 5)

    def test_pipeline(self):
        run_pipeline()
        self.check_files()
