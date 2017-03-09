# Script to launch diffbind package and plot all the figures
# author: oleg.shpynov@jetbrains.com

# Function to load or install library from bioconductor
load_library <- function(name) {
    if (!require(name)) {
        source("http://bioconductor.org/biocLite.R")
        biocLite(name)
    }
    library(name)
}

load_library(DiffBind)
load_library(ggplot2)
load_library(stringr)

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
