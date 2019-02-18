import gzip


def read_fastq_entry(file):
    res = []
    for i in range(0, 4):
        res.append(file.readline())
    return res


def write_fastq_entry(file, entry):
    for line in entry:
        file.write(line)


def g_content(fastq_entry):
    sequence = fastq_entry[1]
    count = 0
    for c in sequence:
        if c == "G":
            count += 1
    # sequence strings are reported with trailing \n, hence -1
    return count * 1.0 / (len(sequence) - 1)


def filter_g_content(input1, input2, output1, output2, threshold=0.9):
    with gzip.open(input1, 'r') as i1:
        with gzip.open(input2, 'r') as i2:
            with gzip.open(output1, 'w') as o1:
                with gzip.open(output2, 'w') as o2:
                    while True:
                        entry1 = read_fastq_entry(i1)
                        entry2 = read_fastq_entry(i2)
                        if entry1[1] == "" or entry2[1] == "":
                            break
                        if g_content(entry1) < threshold and g_content(entry2) < threshold:
                            write_fastq_entry(o1, entry1)
                            write_fastq_entry(o2, entry2)


rule filter_failed_fastq:
    input:
        i1="{anywhere}/{sample}_1.fastq.gz",
        i2="{anywhere}/{sample}_2.fastq.gz"
    output:
        o1="{anywhere}/filtered/{sample}_1.fastq.gz",
        o2="{anywhere}/filtered/{sample}_2.fastq.gz"
    run: filter_g_content(input.i1, input.i2, output.o1, output.o2)


rule filter_failed_fastq_fastqc:
    input: "trimmed/filtered/{sample}.fastq.gz"
    output:
        html="qc/fastqc_filtered/{sample}_fastqc.html",
        zip="qc/fastqc_filtered/{sample}_fastqc.zip"
    log: "logs/fastqc_filtered/{sample}.log"
    wrapper: "0.31.1/bio/fastqc"


rule filter_failed_fastq_multiqc:
    input: expand("qc/fastqc_filtered/{sample}_fastqc.zip", sample=fastq_names())
    output: "multiqc/fastqc_filtered/multiqc.html"
    log: "multiqc/fastqc_filtered/multiqc.log"
    wrapper: "0.31.1/bio/multiqc"


rule filter_failed_fastq_bowtie2:
    input:
        sample=["trimmed/filtered/{sample}_1.fastq.gz", "trimmed/filtered/{sample}_2.fastq.gz"]
    output: "mapped/filtered/{sample}.bam"
    log: "logs/bowtie2_filtered/{sample}.log"
    params:
        index=config["indexes"],
        extra="-X 2000 --dovetail"
    threads: 8
    wrapper: "0.31.1/bio/bowtie2/align"


rule filter_failed_fastq_multiqc_bowtie2:
    input: expand("logs/bowtie2_filtered/{sample}.log", sample=fastq_aligned_names())
    output: "multiqc/bowtie2_filtered/multiqc.html"
    log: "multiqc/bowtie2_filtered/multiqc.log"
    wrapper: "0.31.1/bio/multiqc"