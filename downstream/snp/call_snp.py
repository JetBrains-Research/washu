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
    tmp_dir = pooled_dir / "tmp"
    if not tmp_dir.exists():
        tmp_dir.mkdir()
    bam_name = bam_path.name
    print(bam_path)

    cp = subprocess.run(["samtools view {} | head -n1".format(bam_path)], stdout=subprocess.PIPE, shell=True)
    result = cp.stdout.decode("utf-8")
    parts = result.split('\t')[0].split(':')

    g_id = parts[0]
    g_pu = parts[2]
    g_lb = parts[3]

    fname = os.path.splitext(bam_name)[0] + "_fixed.bam"
    result_bam = tmp_dir / fname
    result_bai = tmp_dir / (fname + ".bai")

    script = "/mnt/stripe/washu/downstream/snp/preprocess_bam.sh"
    cmd = ["bash", script, str(bam_path), str(result_bam), g_id, g_pu, g_lb, donor]
    subprocess.check_call(cmd)

    return result_bam, result_bai


def call_donor_snp(snp_path, ucsc_path, dbsnp_path, donor):
    pooled_dir = snp_path / "pooled"

    if not pooled_dir.exists():
        pooled_dir.mkdir()

    result_path = pooled_dir / "{}.g.vcf.gz".format(donor)
    result_idx_path = pooled_dir / "{}.g.vcf.gz.tbi".format(donor)

    if result_path.exists():
        return result_path

    all_bams = []

    for s in subdirs:
        dir = chip_seq_path / s
        all_bams.extend(dir.glob(donor + "_*_unique.bam"))

    fixed_bams = []

    for bam in all_bams:
        fixed_bams.append(preprocess_bam(pooled_dir, bam, donor))

    cmd = ["/mnt/stripe/tools/gatk-4.0.3.0/gatk", "--java-options", "-Xmx32g", "HaplotypeCaller"]
    cmd += ["-R", str(ucsc_path)]
    cmd += ["--dbsnp", str(dbsnp_path)]

    for bam, _ in fixed_bams:
        cmd += ["-I", str(bam)]

    tmp_path = pooled_dir / "{}_tmp.g.vcf.gz".format(donor)
    tmp_idx_path = pooled_dir / "{}_tmp.g.vcf.gz.tbi".format(donor)

    cmd += ["-O", str(tmp_path), "-ERC", "GVCF"]

    subprocess.check_call(cmd)

    for bam, bai in fixed_bams:
        bam.unlink()
        bai.unlink()

    tmp_path.rename(result_path)
    tmp_idx_path.rename(result_idx_path)

    return result_path
