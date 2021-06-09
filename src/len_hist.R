library(ggplot2)
library(dplyr)
# library(tidyr)   # replace_na
# library(tibble)  # column_to_rownames

###

# NAME <- 'H3K4me3_A549.ENCFF573MUH.hg19'
#NAME <- 'H3K4me3_A549.ENCFF573MUH.hg38'
NAME <- 'H3K4me3_A549.ENCFF832EOL.hg19'
# NAME <- 'H3K4me3_A549.ENCFF832EOL.hg38'
OUT_DIR <- 'Results/'

###

bed_df <- read.delim(paste0('data/', NAME, '.bed'), as.is = TRUE, header = FALSE)
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
