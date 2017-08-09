require("ggplot2")

paths <- commandArgs(TRUE)

read_all <- function(file_name, header = T) {
  merged = NULL

  for (folder in paths) {
    full_path = paste(folder, "/report/", file_name, sep="")
    dt = read.table(file = full_path, sep=',', header = header)
    if (is.null(merged)) {
      merged = cbind(tool = basename(folder), dt)
    } else {
      merged = rbind(merged, cbind(tool = basename(folder), dt))
    }
  }

  return(merged)
}


merged <- read_all("peaks_length.csv")

my_plot <- ggplot(merged, aes(track, length, colour=tool)) + geom_boxplot() + scale_y_log10() + facet_grid(. ~ tool) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave("peaks_length.png", plot=my_plot, width = 12, height = 7)

merged <- read_all("frip_table.csv")

my_plot <- ggplot(merged, aes(tool, peaks)) + geom_jitter(aes(colour=tool)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave("peaks_number.png", plot=my_plot)

my_plot <- ggplot(merged, aes(tool, frip)) + geom_jitter(aes(colour=tool)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave("frip.png", plot=my_plot)
