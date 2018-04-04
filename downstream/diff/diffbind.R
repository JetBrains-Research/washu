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

require_or_install("DiffBind", bioc = TRUE)
require_or_install("ggplot2")
require_or_install("stringr")

# SeeDiffBind DiffBind-globals.Rd for the full list
#
# DBA_SCORE_READS
# dba.count score is number of reads in ChIP
#
# DBA_SCORE_READS_FOLD
# dba.count score is number of reads in ChIP divided by number of reads in Control
#
# DBA_SCORE_READS_MINUS
# dba.count score is number of reads in ChIP minus  number of reads in Control
#
# DBA_SCORE_RPKM
# dba.count score is RPKM of ChIP 
#
# DBA_SCORE_RPKM_FOLD
# dba.count score is RPKM of ChIP divided by RPKM of Control
#
# DBA_SCORE_TMM_READS_FULL
# dba.count score is TMM normalized (using edgeR), using ChIP read counts and Full Library size
#
# DBA_SCORE_TMM_READS_EFFECTIVE
# dba.count score is TMM normalized (using edgeR), using ChIP read counts and Effective Library size
#
# DBA_SCORE_TMM_MINUS_FULL
# dba.count score is TMM normalized (using edgeR), using ChIP read counts minus Control read counts and Full Library size
#
# DBA_SCORE_TMM_MINUS_EFFECTIVE
# dba.count score is TMM normalized (using edgeR), using ChIP read counts minus Control read counts and Effective Library size
#
# DBA_SCORE_TMM_READS_FULL_CPM
# dba.count score is TMM normalized (using edgeR), using ChIP read counts and Full Library size, reported in counts-per-million.
#
# DBA_SCORE_TMM_READS_EFFECTIVE_CPM
# dba.count score is TMM normalized (using edgeR), using ChIP read counts and Effective Library size, reported in counts-per-million.
#
# DBA_SCORE_TMM_MINUS_FULL_CPM
# dba.count score is TMM normalized (using edgeR), using ChIP read counts minus Control read counts and Full Library size, reported in counts-per-million.
#
# DBA_SCORE_TMM_MINUS_EFFECTIVE_CPM
# dba.count score is TMM normalized (using edgeR), using ChIP read counts minus Control read counts and Effective Library size, reported in counts-per-million.
# SCORE_FUNCTIONS = c(DBA_SCORE_READS,
#                     DBA_SCORE_READS_MINUS,
#                     DBA_SCORE_READS_FOLD,
#                     DBA_SCORE_RPKM,
#                     DBA_SCORE_RPKM_FOLD,
#                     DBA_SCORE_TMM_READS_FULL,
#                     DBA_SCORE_TMM_READS_EFFECTIVE,
#                     DBA_SCORE_TMM_MINUS_FULL,
#                     DBA_SCORE_TMM_MINUS_EFFECTIVE,
#                     DBA_SCORE_TMM_READS_FULL_CPM,
#                     DBA_SCORE_TMM_READS_EFFECTIVE_CPM,
#                     DBA_SCORE_TMM_MINUS_FULL_CPM,
#                     DBA_SCORE_TMM_MINUS_EFFECTIVE_CPM)

# SCORE_FUNCTIONS_NAMES =   c('DBA_SCORE_READS',
#                             'DBA_SCORE_READS_MINUS',
#                             'DBA_SCORE_READS_FOLD',
#                             'DBA_SCORE_RPKM',
#                             'DBA_SCORE_RPKM_FOLD',
#                             'DBA_SCORE_TMM_READS_FULL',
#                             'DBA_SCORE_TMM_READS_EFFECTIVE',
#                             'DBA_SCORE_TMM_MINUS_FULL',
#                             'DBA_SCORE_TMM_MINUS_EFFECTIVE',
#                             'DBA_SCORE_TMM_READS_FULL_CPM',
#                             'DBA_SCORE_TMM_READS_EFFECTIVE_CPM',
#                             'DBA_SCORE_TMM_MINUS_FULL_CPM',
#                             'DBA_SCORE_TMM_MINUS_EFFECTIVE_CPM')

SCORE_FUNCTIONS = c(DBA_SCORE_TMM_MINUS_FULL)
SCORE_FUNCTIONS_NAMES =   c('DBA_SCORE_TMM_MINUS_FULL')

# REMOVE_DUPLICATES = c(TRUE, FALSE)
# INSERT_SIZES = c(125, 150) # Default insertSize

REMOVE_DUPLICATES = c(TRUE)
INSERT_SIZES = c(125) # Default insertSize

main <- function(path, insertSize) {
    if (! missing(insertSize)) {
        INSERT_SIZES = c(insertSize)
    }
    print(paste("DiffBind version", packageVersion("DiffBind")))
    print(paste("Processing file", path))

    yo = dba(sampleSheet = path)
    for (insertSize in INSERT_SIZES) {
        for (removeDuplicates in REMOVE_DUPLICATES) {
            for (i in 1 : length(SCORE_FUNCTIONS)) {
                scoreFunction = SCORE_FUNCTIONS[i]
                scoreFunctionName = SCORE_FUNCTIONS_NAMES[i]
                print(paste('INSERT_SIZE', insertSize))
                print(paste('REMOVE_DUPLICATES', removeDuplicates))
                print(paste('SCORE_FUNCTION', scoreFunctionName))
                id = paste('dedup', removeDuplicates, 'f', insertSize, scoreFunctionName, sep = '_')

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

                print(paste('PCA, MA, Volcano plots', id))
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

                # dba.plotVolcano fails if not difference found
                # vlc_pdf = str_replace(path, '.csv', paste('_', id, '_vlc.pdf', sep = ''))
                # pdf(file = vlc_pdf)
                # dba.plotVolcano(yo_contrast)
                # dev.off()
                # print(vlc_pdf)

                print(paste('Save difference to resulting csv file', id))
                db = dba.report(yo_contrast)
                difference_csv = str_replace(path, '.csv', paste('_', id, '_difference.csv', sep = ''))
                write.table(db, difference_csv, sep = ",", row.names = FALSE)
                print(difference_csv)
            }
        }
    }
}

if (! interactive()) {
    args <- commandArgs(TRUE)
    if (length(args) < 1) {
        print("Usage: [executable] config.csv [insertSize]", stderr())
        q(status = 1)
    } else {
        do.call(main, as.list(args))
        warnings()
    }
}
