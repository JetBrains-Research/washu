from pathlib import Path

import subprocess

import os

chip_seq_path = Path("/mnt/stripe/bio/experiments/aging/chipseq")

subdirs = [
    "k27ac/k27ac_20vs20_bams",
    "k27ac/k27ac_20vs20_redo_bams",
    "k27me3/k27me3_20vs20_bams",
    "k36me3/k36me3_20vs20_bams",
    "k4me1/k4me1_20vs20_reseq_bams",
    "k4me3/k4me3_20vs20_bams"]


def preprocess_bam(pooled_dir, bam_path, donor):
    bam_name = bam_path.name
    print(bam_path)

    cp = subprocess.run(["samtools view {} | head -n1".format(bam_path)], stdout=subprocess.PIPE, shell=True)
    result = cp.stdout.decode("utf-8")
    parts = result.split('\t')[0].split(':')

    g_id = parts[0]
    g_pu = parts[2]
    g_lb = parts[3]

    fname = os.path.splitext(bam_name)[0] + "_fixed.bam"
    result_bam = pooled_dir / fname

    script = "/mnt/stripe/washu/downstream/snp/preprocess_bam.sh"
    cmd = ["bash", script, str(bam_path), str(result_bam), g_id, g_pu, g_lb, donor]
    subprocess.check_call(cmd, )


def make_preprocessed_bams_for_donor(snp_path, donor):
    pooled_dir = snp_path / "pooled"

    if not pooled_dir.exists():
        pooled_dir.mkdir()

    all_bams = []

    for s in subdirs:
        dir = chip_seq_path / s
        all_bams.extend(dir.glob(donor + "_*_unique.bam"))

    for bam in all_bams:
        preprocess_bam(pooled_dir, bam, donor)
