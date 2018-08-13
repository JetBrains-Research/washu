library(readr)
library(dplyr)
library(stringr)
library(argparse)

# parse command line arguments
parser <- ArgumentParser()
parser$add_argument('-d', '--inputDirectory', help = 'input directory', required = TRUE)
parser$add_argument('-p', '--prefix', help = 'file name prefix', required = TRUE)
parser$add_argument('-o', '--output', help = 'output file name', required = TRUE)
args <- parser$parse_args()

# get file info
input_directory <- args$inputDirectory
output_file <- args$output
prefix <- args$prefix

meth_files <- list.files(input_directory)
is_initializes <- FALSE
for (file in meth_files) {
  print(file)
  tag <- str_split(file, '\\.', simplify = T)[1, 3]
  if (is_initializes) {
    current_df <- read_tsv(paste0(input_directory, '/', file)) %>% 
      select(chrBase, coverage)
    tagged_names <- c('chrBase', paste0('coverage.', tag))
    colnames(current_df) <- tagged_names
    df <- full_join(df, current_df, by = 'chrBase')
  } else {
    df <- read_tsv(paste0(input_directory, '/', file)) %>% 
      select(chrBase, coverage)
    tagged_names <- c('chrBase', paste0('coverage.', tag))
    colnames(df) <- tagged_names
    is_initializes <- TRUE
  }
}

write_tsv(df, output_file)
