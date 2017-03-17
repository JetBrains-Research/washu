# Script to launch diffbind package and plot all the figures
# author: oleg.shpynov@jetbrains.com

# Function to load or install library from bioconductor
installed <- rownames(installed.packages())

require_or_install <- function(..., bioc = FALSE) {
  for (p in list(...)) {
    if (!(p %in% installed)) {
      if (bioc) {
        if (!exists("biocLite")) {
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

main <- function(path) {
  write(paste("Processing file", path))
  result_pdf = str_replace(path, ".csv", "_result.pdf")
  pdf(file=result_pdf)
  yo=dba(sampleSheet=path)
  yo = dba.count(yo)
  # Plot histogram
  plot(yo)
  # Analyze contrast factor and compute difference
  yo = dba.contrast(yo)
  yo = dba.analyze(yo)
  # Print number of overlap peaks by donors
  olap.rate=dba.overlap(yo,mode=DBA_OLAP_RATE)
  plot(olap.rate,type="b", ylab="# peaks", xlab="Overlap at least this many peaksets")
  # PCA, MA plots
  dba.plotPCA(yo)
  dba.plotMA(yo)
  # Save difference to resulting csv file
  db = dba.report(yo)
  result_csv = str_replace(path, ".csv", "_result.csv")
  write.table(db, result_csv, sep = ",", row.names = FALSE)
  write(paste("Saved", result_csv))
  # IMPORTANT: plot difference histogram on the last step as it fails if there is no difference
  dba.plotBox(yo)
  dev.off()
  write(paste("Saved plots", result_pdf))
}

if (!interactive()) {
  args <- commandArgs(TRUE)
  if (length(args) != 1) {
    write("Usage: [executable] diffbind.csv", stderr())
    q(status = 1)
  } else {
    do.call(main, as.list(args))
    warnings()
  }
}
