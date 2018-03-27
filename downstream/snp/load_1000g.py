import os.path
import subprocess

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
    download_file(base_url + "/" + phenotype, os.path.join(path, phenotype))
    for file in file_names:
        download_file(base_url + file, os.path.join(path, file))
        download_file(base_url + file + ".tbi", os.path.join(path, file + ".tbi"))


def get_1000g(base_path):
    path = os.path.join(base_path, "1000g")
    if not os.path.exists(path):
        os.mkdir(path)

    os.chdir(path)

    result_path = "{}/g1000.vcf.gz".format(path)

    if os.path.exists(result_path):
        return

    download_all(path)


    files = [os.path.join(path, f) for f in file_names]

    tmp_path = ("{}/g1000_tmp.vcf.gz".format(path))
    subprocess.check_call(["bcftools", "concat",
                           "-o", tmp_path,
                           "-O", "z"] + files)


    os.rename(tmp_path, result_path)

    subprocess.check_call(["tabix", "-p", "vcf", result_path])

    for file in file_names:
        os.remove(os.path.join(path, file))
        os.remove(os.path.join(path, file + ".tbi"))