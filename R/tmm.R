# Process reads normalization using TMM from edgeR package
# author: oleg.shpynov@jetbrains.com

# Function to load or install library from bioconductor
installed <- rownames(installed.packages())

require_or_install <- function(..., bioc = FALSE) {
    for (p in list(...)) {
        if (! (p %in% installed)) {
            if (bioc) {
                if (! exists("biocLite")) {
                    source("http://bioconductor.org/biocLite.R")
                }
                biocLite(p)
            } else {
                install.packages(p, repos = "http://cran.us.r-project.org")
            }
        }

        require(p, character.only = TRUE)
    }
}
require_or_install("edgeR")

normalize <- function(input_path, libsizes_path, output_path) {
    df = read.table(input_path, header = TRUE)
    libsizes = read.table(libsizes_path)
    names = colnames(df)
    groups = rep('', ncol(df))
    for (i in 1:ncol(df)) { if (startsWith(names[i], "O")) groups[i] = 'O' else groups[i] = 'Y' }
    # Transpose libsizes, so that we can do indexing
    libsizes.T = t(libsizes[,2])
    colnames(libsizes.T) <- libsizes[,1]

    # Prepare list of libsizes
    ls=rep(0, ncol(df))
    for (i in 1:ncol(df)) {ls[i] = libsizes.T[1,names[i]]}

    # Create EdgeR object and use TMM normalization
    # See analyse.R:1380
    res = DGEList(df, lib.size = ls, group = groups)
    res = calcNormFactors(res,method="TMM",doWeighting=F)

    # Update counts given normalization factors
    counts <- res$counts
    sizes  <- res$samples$lib.size * res$samples$norm.factors
    counts <- t(t(counts)/sizes)
    write.table(counts, output_path, row.names = F, quote = F, sep='\t')
}


if (!interactive()) {
    args <- commandArgs(TRUE)
    if (length(args) != 3) {
        print("Usage: [executable] to_normalize.tsv libraries.tsv normalized.tsv", stderr())
        q(status = 1)
    } else {
        do.call(normalize, as.list(args))
        warnings()
    }
}