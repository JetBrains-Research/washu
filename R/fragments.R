# Script to plot insert size metrics for alignment.
library(ggplot2)

main <- function(insertSizeMetrics, fragments) {
  inserts <- scan(insertSizeMetrics)
  inserts <- table(inserts)
  inserts <- as.data.frame(inserts)

  ggplot(data=inserts, aes(x=inserts, y=Freq, group=1)) +
    geom_line() + theme_bw() + scale_x_discrete(breaks=c(0, 100, 200, 300, 400, 500)) +
    xlab("fragments length") + ylab("frequency")
  ggsave(fragments, width=10, height=5, units="in")
  print(sprintf("Processed %s, result: %s", insertSizeMetrics, fragments))
}

if (!interactive()) {
  args <- commandArgs(TRUE)
  if (length(args) != 2) {
    write("Usage: [executable] <InsertSizeMetrics.txt> <Fragments.png>",
          stderr())
    q(status = 1)
  } else {
    main(args[1], args[2])
    warnings()
  }
}