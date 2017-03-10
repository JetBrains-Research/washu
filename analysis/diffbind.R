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
  yo=dba(sampleSheet=path)
  yo = dba.count(yo)
  plot(yo)
  yo = dba.contrast(yo)
  yo = dba.analyze(yo)
  dba.plotMA(yo)
  dba.plotPCA(yo,contrast=1)
  pvals=dba.plotBox(yo)
  olap.rate=dba.overlap(yo,mode=DBA_OLAP_RATE)
  plot(olap.rate,type="b", ylab="# peaks", xlab="Overlap at least this many peaksets")
  db = dba.report(yo)
  result = str_replace(path, ".csv", "_result.csv")
  write(paste("Saved", result))
  write.table(db, result, sep = ",", row.names = FALSE)
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
