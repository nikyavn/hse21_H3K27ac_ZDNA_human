library(ggplot2)
library(dplyr)
# library(tidyr)   # replace_na
# library(tibble)  # column_to_rownames

###

# NAME <- 'H3K27ac_H9.ENCFF997MGG.hg19'
#NAME <- 'H3K27ac_H9.ENCFF997MGG.hg38'
#NAME <- 'H3K27ac_H9.ENCFF365GJO.hg19'
 NAME <- 'H3K27ac_H9.ENCFF365GJO.hg38'
OUT_DIR <- 'D:/HSE/3_course/bioinf/project/Results'

###

bed_df <- read.delim(paste0('D:/HSE/3_course/bioinf/project/bed_files/', NAME, '.bed'), as.is = TRUE, header = FALSE)
colnames(bed_df) <- c('chrom', 'start', 'end', 'name', 'score')
bed_df$len <- bed_df$end - bed_df$start
head(bed_df)

# hist(bed_df$len)

ggplot(bed_df) +
  aes(x = len) +
  geom_histogram() +
  ggtitle(NAME, subtitle = sprintf('Number of peaks = %s', nrow(bed_df))) +
  theme_bw()
ggsave(paste0('len_hist.', NAME, '.pdf'), path = OUT_DIR)
