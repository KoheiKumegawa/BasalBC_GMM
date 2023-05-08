#------------------------------------------------------------------------------
# RNA_preprocess.R
#------------------------------------------------------------------------------
#Kohei Kumegawa, Japanese Foundation for Cancer Research
library(SummarizedExperiment)
library(data.table)
library(dplyr)
library(ggplot2)
countToTpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

#make summrized experiment
counts <- fread("data/counts_exon_uns.txt")
len <- counts$Length
gene <- counts$Geneid
counts <- counts[, -c(1:6)] %>% 
  `colnames<-`(., gsub("_hg38_Aligned.sortedByCoord.out.bam", "", colnames(counts[, -c(1:6)]))) %>%
  as.matrix %>% `rownames<-`(., gene) 
tpms <- countToTpm(counts, len)
se <- SummarizedExperiment(assays = list(counts = counts, log2tpm = log2(tpms+1)),
                           rowData = DataFrame(symbol = gene),
                           colData = DataFrame(sample = colnames(counts)))

colData(se)$cell <- stringr::str_split(colData(se)$sample, "-", simplify = T)[,2]
colData(se)$siType <- stringr::str_split(colData(se)$sample, "-", simplify = T)[,3]
colData(se)$siRep <- stringr::str_split(colData(se)$sample, "-", simplify = T)[,4]
colData(se)$siRep[c(5,6,11,12)] <- 1
colData(se)$expRep <- rep(c(1,2),6)
colData(se)$sampleType <- paste0(colData(se)$cell, "_", colData(se)$siType, "_", colData(se)$siRep)

saveRDS(se, "rds/se.rds")
