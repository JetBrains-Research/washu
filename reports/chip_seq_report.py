#!/usr/bin/env python

import os
import sys
import pandas as pd
import subprocess

help_data = """
Usage: chip_seq_report.py [chrom.sizes] [peaks files]

Script creates report folder with ChIP-seq peaks statistics:    
 counts_stat.csv - consensus table. Table contains statistics for coverage of all tracks peaks.   
 frip_table.csv  - table with FRiPs 
 length.png      - distribution of peaks length
 jaccard.csv     - jaccard statistics for all pairs of peaks

"""


def write_frips():
    files = [file for file in os.listdir(".") if file.endswith("_rip.csv")]
    with open(os.path.join("report", "frip_table.csv"), "w") as result:
        result.write("file,peaks,frip\n")
        for file in files:
            data = pd.read_csv(file)
            file = data["file"].loc[0]
            reads = float(data["reads"].loc[0])
            peaks = int(data["peaks"].loc[0])
            rip = float(data["rip"].loc[0])
            result.write("{},{},{}\n".format(file, peaks, rip * 100 / reads))


def process(chrom_sizes, peaks):
    if not os.path.exists("report"):
        os.mkdir("report")

    print("writing FRiP table")
    write_frips()

    print("writing peaks length")
    write_peak_length(peaks)

    print("writing intersection counts")
    count_intersection(chrom_sizes, peaks)

    print("writing jaccard matrix")
    with open(os.path.join("report", "jaccard.csv"), "w") as result:
        result.write("f1,f2,jaccard\n")
        for p1 in peaks:
            for p2 in peaks:
                output = subprocess.check_output("bedtools jaccard -a {} -b {}".format(p1, p2), shell=True)
                line = output.split("\n")[1]
                result.write("{},{},{}\n".format(p1, p2, float(line.split("\t")[2])))

    make_plots_script()

    print("done")


def count_intersection(chrom_sizes, peaks):
    command = "cat {} | sort -k 1,1 | bedtools genomecov -bg -i - -g {} >{}".format(
        " ".join(peaks), chrom_sizes, os.path.join("report", "counts.bed"))
    os.system(command)
    counts = [0] * len(peaks)
    for line in read_all_lines(os.path.join("report", "counts.bed")):
        parts = line.split("\t")
        count = int(parts[3])
        length = int(parts[2]) - int(parts[1])

        counts[count - 1] += length
    s = sum(counts)
    with open(os.path.join("report", "counts_stat.csv"), "w") as result:
        for i in range(len(counts)):
            result.write("{},{},{}\n".format(i + 1, float(counts[i]) / s, counts[i]))


def write_peak_length(peaks):
    with open(os.path.join("report", "peaks_length.csv"), "w") as result:
        result.write("track,length\n")
        for i, peak_file in enumerate(peaks):
            for line in read_all_lines(peak_file):
                parts = line.split("\t")
                length = int(parts[2]) - int(parts[1])
                result.write("{},{}\n".format(i, length))


r_file_text = """

require("ggplot2")

dt = read.csv("peaks_length.csv")
plt = ggplot(dt, aes(length)) + scale_x_log10() + geom_histogram(bins=40) + facet_grid(track ~ .)
ggsave(filename="length.png", plot=plt)

consensus.df = read.table("./counts.bed", header = FALSE)
consensus.scores = sapply(strsplit(as.character(consensus.df$V4), split = "\\|"), length)
consensus = min(consensus.scores):max(consensus.scores)
cons.numbers = sapply(consensus, function(x) sum(consensus.scores == x))
plt2 = ggplot(data.frame(consensus.peaks = cons.numbers), aes(x="", y=consensus.peaks, fill=consensus)) + 
    geom_bar(width=1, stat = "Identity") + coord_polar("y", start = 0) + ggtitle("Consensus by peaks")
ggsave(filename = "./consensus.png", plot = plt2)

"""


def make_plots_script():
    script_name = "make_plots.r"
    with open(os.path.join("report", script_name), "w") as plots:
        plots.write(r_file_text)
    p = subprocess.Popen(["Rscript", script_name], cwd="report")
    p.wait()


def read_all_lines(peak_file):
    with open(peak_file, "r") as f:
        return f.readlines()


def main():
    args = sys.argv

    if len(args) < 2:
        print(help_data)
        sys.exit(1)
    process(args[1], args[2:])


if __name__ == "__main__":
    main()
