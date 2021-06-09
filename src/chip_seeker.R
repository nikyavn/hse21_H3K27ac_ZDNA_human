
 #Installing BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ChIPseeker")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene", force=TRUE)
BiocManager::install("clusterProfiler", force=TRUE)
BiocManager::install("GenomicFeatures", force=TRUE)
BiocManager::install("org.Hs.eg.db", force=TRUE)
BiocManager::install("ChIPpeakAnno", force=TRUE)

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
###

NAME <- 'H3K27ac_H9.ENCFF365GJO.hg19.filtered'
#NAME <- 'H3K27ac_H9.ENCFF997MGG.hg19.filtered'
DATA_DIR<- 'D:/HSE/3_course/bioinf/project/bed_files/'
BED_FN <- paste0(DATA_DIR, NAME, '.bed')
OUT_DIR <- 'D:/HSE/3_course/bioinf/project/Results/'

###

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

peakAnno <- annotatePeak(BED_FN, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

png(paste0(OUT_DIR, 'chip_seeker.', NAME, '.plotAnnoPie.png'))
plotAnnoPie(peakAnno)
dev.off()
