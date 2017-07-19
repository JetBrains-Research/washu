# Script to analysize MACS2 in terms of peaks characteristics.
#
# Use the following BASH script to create report.csv file required as input.
#
# head -1 peaks_k4me1_macs_broad_0.1/macs2_report.csv > report.csv
# find . -name '*report.csv' -depth 2 | xargs cat | grep -v 'sample,' | grep -v 'unique' | awk -v FS=',' '{ printf("%s,%s_NE,PS_%s_PE\n",$0,$1,$1) }' |\
# sed -e 's#_hg19[^,]*_NE##g' -e 's#PS_[^,]*hg19_##g' -e 's#_macs2[^,]*_PE##g' >> report.csv

library(ggplot2)
library(stringr)
library(dplyr)

plot_distribution = function(df_lens, title) {
  plot(ggplot(df_lens, aes(x=len, colour=params)) + geom_density() + ggtitle(title))
  # Hide x tick marks, labels, and grid lines
  plot(ggplot(df_lens, aes(x=params, y=len, colour=params)) + geom_boxplot() + ggtitle(title) + scale_x_discrete(breaks=NULL))
}

process = function(df) {
  # Summary statistics
  plot(ggplot(df, aes(fill=params)) + geom_bar(aes(x=name, y=frip), stat="identity", position="dodge"))
  plot(ggplot(df, aes(fill=params)) + geom_bar(aes(x=name, y=peaks), stat="identity", position="dodge"))
  plot(ggplot(df, aes(x=peaks, y=frip, shape=name, color=params, size=5)) + geom_point())
  
  
  unique_names=unique(df['name'])
  # TODO(shpynov): is there a better way of iteration over rows?
  for (i in 1:dim(unique_names)[1]) { 
    name = unique_names$name[i]
    print(paste("Processing", name))
    df_name = df[df$name == name, c('sample', 'params')]
    df_lens = NULL
    for (r in 1:dim(df_name)[1]) {
      sample=df_name$sample[r]
      params=df_name$params[r]
      folder = paste('/Users/oleg/Desktop/peaks_k4me1/peaks_k4me1_macs_', params, sep = '')
      file=str_replace(paste(folder, sample, sep='/'), '_macs2.log', '_peaks.broadPeak')
      print(file)
      bed = read.table(file, header = FALSE, sep='\t')
      bed$len = bed[, 3] - bed[, 2]
      bed$params = params
      
      bed=bed[, c('len', 'params')]
      if (is.null(df_lens)) {
        df_lens = bed
      } else {
        df_lens = rbind(df_lens, bed)
      }
    }
    plot_distribution(df_lens, name)
    plot_distribution(df_lens[df_lens$len<10000, ], paste(name, "< 10Kbp"))
    plot_distribution(df_lens[df_lens$len<1000, ], paste(name, "< 1Kbp"))
    plot_distribution(df_lens[df_lens$len<500, ], paste(name, "< 500bp"))
  }
}

main = function() {
  df = read.csv('/Users/oleg/Desktop/peaks_k4me1/report.csv')

  output = '/Users/oleg/Desktop/peaks_k4me1/plots_broad_01.pdf'
  print(paste("Processing as is", output))
  pdf(output)
  process(filter(df, params == "broad_0.1" | 
                   params == "broad_0.1_nolambda" | 
                   params == "broad_0.1_bw150" | 
                   params == "broad_0.1_bw600" |
                   params == "broad_0.1_mfold10-30" |
                   params == "broad_0.1_mfold2-100" |
                   params == "q0.05_broad_0.5"))
  dev.off()
  
  output = '/Users/oleg/Desktop/peaks_k4me1/plots_filtered.pdf'
  print(paste("Processing as is", output))
  pdf(output)
  process(filter(df, params == "broad_0.1" | 
                   params == "broad_0.1_nolambda" |
                   params == "q0.01_broad" |
                   params == "q0.01_broad_0.5" | 
                   params == "q0.05_broad_0.5"))
  dev.off()

  output = '/Users/oleg/Desktop/peaks_k4me1/plots.pdf'
  print(paste("Processing as is", output))
  pdf(output)
  process(df)
  dev.off()
  
}
# Just do it!
main()