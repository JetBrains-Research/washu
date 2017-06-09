# Script to launch diffbind package and plot all the figures
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

require_or_install("DiffBind")
require_or_install("ggplot2")
require_or_install("stringr")

main <- function(path, fragmentSize) {
    print(paste("DiffBind version", packageVersion("DiffBind")))
    print(paste("Processing file", path))
    print(paste("Fragment size", fragmentSize))

    yo = dba(sampleSheet = path)
    yo = dba.count(yo, bRemoveDuplicates=TRUE, fragmentSize=fragmentSize)

    # Write counts table
    counts <- dba.peakset(yo, bRetrieve=TRUE)
    counts_csv = str_replace(path, ".csv", "_counts.csv")
    write.table(counts, counts_csv, sep = ",", row.names = FALSE)

    # Plot histogram
    pdf(file = str_replace(path, ".csv", "_histogram.pdf"))
    plot(yo)
    dev.off()

    # Analyze contrast factor and compute difference
    yo = dba.contrast(yo)
    yo = dba.analyze(yo)

    # Print number of overlap peaks by donors
    pdf(file = str_replace(path, ".csv", "_overlap.pdf"))
    overlap_rate = dba.overlap(yo, mode = DBA_OLAP_RATE)
    plot(overlap_rate, type = "b", ylab = "# peaks", xlab = "Overlap at least this many peaksets")
    dev.off()

    # PCA, MA plots
    pdf(file = str_replace(path, ".csv", "_pca.pdf"))
    dba.plotPCA(yo)
    dev.off()

    pdf(file = str_replace(path, ".csv", "_ma.pdf"))
    dba.plotMA(yo)
    dev.off()

    # Save difference to resulting csv file
    db = dba.report(yo)
    result_csv = str_replace(path, ".csv", "_result.csv")
    write.table(db, result_csv, sep = ",", row.names = FALSE)
    write(paste("Saved", result_csv))

    # IMPORTANT: plot difference on the last step as it fails if there is no difference
    pdf(file = str_replace(path, ".csv", "_difference_histogram.pdf"))
    plot(yo, contrast = 1)
    dev.off()

    pdf(file = str_replace(path, ".csv", "_difference_pca.pdf"))
    dba.plotPCA(yo, contrast = 1)
    dev.off()

    pdf(file = str_replace(path, ".csv", "_difference.pdf"))
    dba.plotBox(yo)
    dev.off()
}

if (! interactive()) {
    args <- commandArgs(TRUE)
    if (length(args) != 2) {
        write("Usage: [executable] diffbind.csv fragment_size", stderr())
        q(status = 1)
    } else {
        do.call(main, as.list(args))
        warnings()
    }
}
