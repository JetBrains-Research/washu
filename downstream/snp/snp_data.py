import os
import subprocess
from pathlib import Path

base_url = "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/"

file_names = [
    "ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
    "ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
    "ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
    "ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
    "ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
    "ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
    "ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
    "ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
    "ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
    "ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
    "ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
    "ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
    "ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
    "ALL.chr14.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
    "ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
    "ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
    "ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
    "ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
    "ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
    "ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
    "ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
    "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
    "ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz"]


def download_file(url, result_file):
    print("Downloading file: {}".format(url))
    subprocess.check_call(["wget", "-c", "-O", result_file, url])


def download_all(path):
    phenotype = "integrated_call_samples.20130502.ALL.ped"
    download_file(base_url + "/" + phenotype, path / phenotype)
    for file in file_names:
        download_file(base_url + file, path / file)
        download_file(base_url + file + ".tbi", path / (file + ".tbi"))


def get_1000g(snp_path):
    path = snp_path / "1000g"
    if not path.exists():
        path.mkdir()

    result_path = path / "g1000.vcf.gz"

    if result_path.exists():
        print("1000 genomes exists: {}".format(result_path))
        return

    path.chdir()

    download_all(path)

    files = [str(path / f) for f in file_names]

    tmp_path = path / "g1000_tmp.vcf.gz"
    subprocess.check_call(["bcftools", "concat",
                           "-o", tmp_path,
                           "-O", "z"] + files)

    tmp_path.rename(result_path)

    subprocess.check_call(["tabix", "-p", "vcf", result_path])

    for file in file_names:
        (path / file).unlink()
        (path / (file + ".tbi")).unlink()


def get_snp_path():
    snp_dir = Path("/mnt/stripe/bio/experiments/snp")
    if not snp_dir.exists():
        snp_dir.mkdir()

    return Path(snp_dir)


def get_ucsc_path():
    snp_path = get_snp_path()
    ucsc_dir = snp_path / "ucsc"
    if ucsc_dir.exists():
        return ucsc_dir

    ucsc_tmp = snp_path / "ucsc_tmp"

    if not ucsc_tmp.exists():
        ucsc_tmp.mkdir()

    return ucsc_dir
