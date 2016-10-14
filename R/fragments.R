library(ggplot2)
inserts <- scan("InsertSizeMetrics.txt")
inserts <- table(inserts)
inserts <- as.data.frame(inserts)

ggplot(data=inserts, aes(x=inserts, y=Freq, group=1)) +
    geom_line() + theme_bw() + scale_x_discrete(breaks=c(0, 100, 200, 300, 400, 500)) +
    xlab("fragments length") + ylab("frequency")
ggsave("frarments.png", width=10, height=5, units="in")