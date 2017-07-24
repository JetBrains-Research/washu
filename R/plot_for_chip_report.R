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

ggplot(merged, aes(track, length, colour=tool)) + geom_boxplot() + scale_y_log10() + facet_grid(. ~ tool)

ggsave("peaks_length.png")

merged <- read_all("frip_table.csv")

ggplot(merged, aes(tool, peaks)) + geom_jitter(aes(colour=tool))

ggsave("peaks_number.png")

ggplot(merged, aes(tool, frip)) + geom_jitter(aes(colour=tool))

ggsave("frip.png")
