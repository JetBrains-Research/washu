from pathlib import Path

import subprocess

import os

chip_seq_path = Path("/mnt/stripe/bio/experiments/aging/chipseq")

gatk_path = "/mnt/stripe/tools/gatk-4.0.3.0/gatk"

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

    cp = subprocess.run(["samtools view {} | head -n1".format(bam_path)], stdout=subprocess.PIPE,
                        shell=True)
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


def call_donor_variants(snp_path, ucsc_path, dbsnp_path, donor):
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

    cmd = [gatk_path, "--java-options", "-Xmx32g", "HaplotypeCaller"]
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


def combine_gvcfs(paths, snp_path, reference_path):
    combined_dir = snp_path / "combined"

    if not combined_dir.exists():
        combined_dir.mkdir()

    result_path = combined_dir / "combined.g.vcf.gz"
    result_idx_path = combined_dir / "combined.g.vcf.gz.tbi"

    if result_path.exists():
        return result_path

    tmp_path = combined_dir / "combined_tmp.g.vcf.gz"
    tmp_idx_path = combined_dir / "combined_tmp.g.vcf.gz.tbi"

    cmd = [gatk_path, "CombineGVCFs", "-R", str(reference_path)]

    for path in paths:
        cmd += ["--variant", str(path)]

    cmd += ["-O", str(tmp_path)]

    print(cmd)

    subprocess.check_call(cmd)

    tmp_path.rename(result_path)
    tmp_idx_path.rename(result_idx_path)


def call_snp(snp_path, reference_path):
    combined_dir = snp_path / "combined"

    if not combined_dir.exists():
        combined_dir.mkdir()

    input_path = combined_dir / "combined.g.vcf.gz"

    result_path = combined_dir / "cohort.g.vcf.gz"
    result_idx_path = combined_dir / "cohort.g.vcf.gz.tbi"

    if result_path.exists():
        return result_path

    tmp_path = combined_dir / "cohort_tmp.g.vcf.gz"
    tmp_idx_path = combined_dir / "cohort_tmp.g.vcf.gz.tbi"

    cmd = [gatk_path, "--java-options", "-Xmx32g",
           "GenotypeGVCFs", "-R", str(reference_path),
           "-V", str(input_path),
           "-O", str(tmp_path)]

    print(cmd)

    subprocess.check_call(cmd)

    tmp_path.rename(result_path)
    tmp_idx_path.rename(result_idx_path)


def filter_snp(snp_path, reference_path):
    combined_dir = snp_path / "combined"
    cohort_path = combined_dir / "cohort.vcf.gz"

    result_path = combined_dir / "final.vcf.gz"
    result_idx_path = combined_dir / "final.vcf.gz.tbi"

    tmp1_path = combined_dir / "output_fg.vcf.gz"
    tmp2_path = combined_dir / "final_tmp.vcf.gz"
    tmp2_idx_path = combined_dir / "final_tmp.vcf.gz"

    cmd = [gatk_path, "VariantFiltration", "-R", str(reference_path),
           "-V", str(cohort_path), "-O", str(tmp1_path),
           "--genotype-filter-expression", "DP < 2", "--genotype-filter-name", "Low_DP",
           "--set-filtered-genotype-to-no-call"]

    subprocess.check_call(cmd)

    cmd = [gatk_path, "SelectVariants", "-R", str(reference_path),
           "--exclude-filtered",
           "--max-nocall-number", "0",
           "-V", str(tmp1_path),
           "-O", str(tmp2_path)]

    subprocess.check_call(cmd)

    tmp2_path.rename(result_path)
    tmp2_idx_path.rename(result_idx_path)
