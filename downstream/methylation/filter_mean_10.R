library(readr)
library(dplyr)
library(argparse)

# parse command line arguments
parser <- ArgumentParser()
parser$add_argument('-i', '--input', help = 'input', required = TRUE)
parser$add_argument('-o', '--output', help = 'output file name', required = TRUE)
args <- parser$parse_args()

input <- args$input
output_file <- args$output

# subset rows: mean(coverage) >= 10
coverage_df <- read_tsv(input)
coverage_matrix <- as.matrix(select(coverage_df, -chrBase))
coverage_matrix[is.na(coverage_matrix)] <- 0
covered_idx <- which(rowMeans(coverage_matrix) >= 10)
coverage_df <- coverage_df[covered_idx, ]
write_tsv(coverage_df, output_file)
