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

SCORE_FUNCTIONS = c(DBA_SCORE_TMM_MINUS_FULL,
                    DBA_SCORE_TMM_MINUS_FULL_CPM,
                    DBA_SCORE_TMM_MINUS_EFFECTIVE,
                    DBA_SCORE_READS_MINUS)
SCORE_FUNCTIONS_NAMES = c('DBA_SCORE_TMM_MINUS_FULL',
                    'DBA_SCORE_TMM_MINUS_FULL_CPM',
                    'DBA_SCORE_TMM_MINUS_EFFECTIVE',
                    'DBA_SCORE_READS_MINUS')
for (i in 1:length(SCORE_FUNCTIONS)) {
    print(paste(SCORE_FUNCTIONS_NAMES[i], SCORE_FUNCTIONS[i]))
}


REMOVE_DUPLICATES = c(TRUE, FALSE)
INSERT_SIZES = c(125) # Default insertSize

main <- function(path, insertSize) {
    if (!missing(insertSize) && ! (insertSize %in% INSERT_SIZES)) {
        INSERT_SIZES = c(insertSize, INSERT_SIZES)
    }
    print(paste("DiffBind version", packageVersion("DiffBind")))
    print(paste("Processing file", path))

    yo = dba(sampleSheet = path)
    for (insertSize in INSERT_SIZES) {
        for (removeDuplicates in REMOVE_DUPLICATES) {
            for (scoreFunction in SCORE_FUNCTIONS) {
                print(paste('INSERT_SIZE', insertSize))
                print(paste('REMOVE_DUPLICATES', removeDuplicates))
                print(paste('SCORE_FUNCTION', scoreFunction))
                id = paste(removeDuplicates, insertSize, scoreFunction, sep = '_')

                print(paste('Count', id))
                yo_counts = dba.count(yo,
                    bRemoveDuplicates = removeDuplicates, fragmentSize = insertSize, score = scoreFunction)
                print(paste('PeakSets', id))
                peaksets <- dba.peakset(yo_counts, bRetrieve = TRUE)
                scores_csv = str_replace(path, '.csv', paste('_', id, '_counts.csv', sep = ''))
                write.table(peaksets, scores_csv, sep = ",", row.names = FALSE)
                print(scores_csv)

                print(paste('Histograms', id))
                hist_pdf = str_replace(path, '.csv', paste('_', id, '_hist.pdf', sep = ''))
                pdf(file = hist_pdf)
                plot(yo_counts)
                dev.off()
                print(hist_pdf)

                print(paste('Analyze contrast factor and compute difference', id))
                yo_contrast = dba.contrast(yo_counts)
                yo_contrast = dba.analyze(yo_contrast)

                print(paste('Number of overlap peaks by donors', id))
                peaks_overlap_pdf = str_replace(path, ".csv", "_overlap.pdf")
                pdf(file = peaks_overlap_pdf)
                overlap_rate = dba.overlap(yo_contrast, mode = DBA_OLAP_RATE)
                plot(overlap_rate, type = "b", ylab = "# peaks", xlab = "Overlap at least this many peaksets")
                dev.off()
                print(peaks_overlap_pdf)

                print(paste('PCA, MA plots', id))
                pca_pdf = str_replace(path, '.csv', paste('_', id, '_pca.pdf', sep = ''))
                pdf(file = pca_pdf)
                dba.plotPCA(yo_contrast)
                dev.off()
                print(pca_pdf)

                ma_pdf = str_replace(path, '.csv', paste('_', id, '_ma.pdf', sep = ''))
                pdf(file = ma_pdf)
                dba.plotMA(yo_contrast)
                dev.off()
                print(ma_pdf)

                print(paste('Save difference to resulting csv file', id))
                db = dba.report(yo_contrast)
                difference_csv = str_replace(path, '.csv', paste('_', id, '_difference.csv', sep = ''))
                write.table(db, difference_csv, sep = ",", row.names = FALSE)
                print(difference_csv)
            }
        }
    }
}

if (!interactive()) {
    args <- commandArgs(TRUE)
    if (length(args) < 1) {
        print("Usage: [executable] config.csv [insertSize]", stderr())
        q(status = 1)
    } else {
        do.call(main, as.list(args))
        warnings()
    }
}
