#------------------------------------------------------------------------------
# ChIP_analysis.R
#------------------------------------------------------------------------------
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ChIPseeker)
source("code/edgeR_PairwiseFunction.R")
'%ni%' <- Negate('%in%')

BT549_se <- readRDS("rds/BT549_K27AC_se_v2.rds")
Hs578T_se <- readRDS("rds/Hs578T_K27AC_se_v2.rds")

#----------------- PCA -----------------# 
pca1 <- prcomp(t(assays(BT549_se)$log2cpm))
df <- data.frame(pca1$x[,c(1,2)], Sample = c("siFOXD1_1", "siFOXD1_1", "siFOXD1_2", "siFOXD1_2", "siNC", "siNC"))
p3 <- ggplot(df, aes(x = PC1, y = PC2, color = Sample)) + geom_point(size = 1) + ArchR::theme_ArchR() + ggtitle("BT549")

pca2 <- prcomp(t(assays(Hs578T_se)$log2cpm))
df <- data.frame(pca2$x[,c(1,2)], Sample = c("siFOXD1_1", "siFOXD1_1", "siFOXD1_2", "siFOXD1_2", "siNC", "siNC"))
p4 <- ggplot(df, aes(x = PC1, y = PC2, color = Sample)) + geom_point(size = 1) + ArchR::theme_ArchR() + ggtitle("Hs578T")

summary(pca1) 
summary(pca2) 

pdf("output/Plots/v2_C_PCA.pdf", width = 3, height = 3.5)
p3
p4
dev.off()

#----------------- Differential analysis -----------------# 
#---- BT549
#siFOXD1_1 vs scNC
idx <- stringr::str_split(colnames(BT549_se), pattern = "-") %>% lapply(., function(i) paste0(i[c(1:3)], collapse = "_")) %>% unlist
colData(BT549_se) <- DataFrame(row.names = colnames(BT549_se), sampleType = idx)

diff_ls <- list(BT549_siFOXD1_1 = list("BT549_siFOXD1_1", "BT549_siNC_1"), BT549_siFOXD1_2 = list("BT549_siFOXD1_2","BT549_siNC_1"))
BT549_DiffTest <- lapply(diff_ls, function(x) edgeR_pairwise(BT549_se, compareCol = "sampleType", topGroup = x[[1]], bottomGroup = x[[2]]))

BT549_df_ls <- lapply(names(BT549_DiffTest), function(x){
  df <- assay(BT549_DiffTest[[x]])[, c("log2FC", "FDR")] %>% data.frame
  df$mlog10FDR <- -log10(df$FDR)
  df$signif <- "N"
  df$signif[which(df$log2FC > 0 & df$FDR < 0.1)] <- "siFOXD1_Up"
  df$signif[which(df$log2FC < 0 & df$FDR < 0.1)] <- "siFOXD1_Dn"
  df$signif <- factor(df$signif, levels = c("N", "siFOXD1_Up", "siFOXD1_Dn"))
  df <- df[order(df$signif),]
  return(df)
})
names(BT549_df_ls) <- names(BT549_DiffTest)

p5 <- lapply(names(BT549_df_ls), function(x){
  df <- BT549_df_ls[[x]]
  out <- ggplot(df, aes(x = log2FC, y = mlog10FDR, color = signif)) + geom_point_rast(size = 0.8) + ArchR::theme_ArchR() +
    geom_hline(yintercept = -log10(0.1), lty = "dashed") + geom_vline(xintercept = c(0), lty = "dashed") + 
    scale_color_manual(values = c("N" = "lightgray", "siFOXD1_Up" = "red", "siFOXD1_Dn" = "blue")) +
    ggtitle(x) + labs(x = "log2FC", y = "-log10(FDR)")
  return(out)
})

pdf("output/Plots/v2_C_BT549_VolcanoPlot.pdf", width = 5, height = 5)
p5
dev.off()

#export
peaks <- rowRanges(BT549_se)
gr1 <- GRangesList(GRangesList(BT549_siFOXD1_1_Up = peaks[rownames(BT549_df_ls$BT549_siFOXD1_1)[BT549_df_ls$BT549_siFOXD1_1$signif == "siFOXD1_Up"]], 
                               BT549_siFOXD1_1_Dn = peaks[rownames(BT549_df_ls$BT549_siFOXD1_1)[BT549_df_ls$BT549_siFOXD1_1$signif == "siFOXD1_Dn"]], 
                               BT549_siFOXD1_2_Up = peaks[rownames(BT549_df_ls$BT549_siFOXD1_2)[BT549_df_ls$BT549_siFOXD1_2$signif == "siFOXD1_Up"]], 
                               BT549_siFOXD1_2_Dn = peaks[rownames(BT549_df_ls$BT549_siFOXD1_2)[BT549_df_ls$BT549_siFOXD1_2$signif == "siFOXD1_Dn"]]))
lapply(names(gr1), function(x){
  g <- gr1[[x]]
  d <- data.frame(seqnames = seqnames(g), start = start(g)-1, end = end(g))
  write.table(d, paste0("output/output_bed_v2//", x, ".bed"), row.names = F, col.names = F, quote = F, sep = "\t")
})

#---- Hs578T
idx <- stringr::str_split(colnames(Hs578T_se), pattern = "-") %>% lapply(., function(i) paste0(i[c(1:3)], collapse = "_")) %>% unlist
colData(Hs578T_se) <- DataFrame(row.names = colnames(Hs578T_se), sampleType = idx)

diff_ls <- list(Hs578T_siFOXD1_1 = list("Hs578T_siFOXD1_1", "Hs578T_siNC_1"), Hs578T_siFOXD1_2 = list("Hs578T_siFOXD1_2","Hs578T_siNC_1"))
Hs578T_DiffTest <- lapply(diff_ls, function(x) edgeR_pairwise(Hs578T_se, compareCol = "sampleType", topGroup = x[[1]], bottomGroup = x[[2]]))

Hs578T_df_ls <- lapply(names(Hs578T_DiffTest), function(x){
  df <- assay(Hs578T_DiffTest[[x]])[, c("log2FC", "FDR")] %>% data.frame
  df$mlog10FDR <- -log10(df$FDR)
  df$signif <- "N"
  df$signif[which(df$log2FC > 0 & df$FDR < 0.1)] <- "siFOXD1_Up"
  df$signif[which(df$log2FC < 0 & df$FDR < 0.1)] <- "siFOXD1_Dn"
  df$signif <- factor(df$signif, levels = c("N", "siFOXD1_Up", "siFOXD1_Dn"))
  df <- df[order(df$signif),]
  return(df)
})
names(Hs578T_df_ls) <- names(Hs578T_DiffTest)

p6 <- lapply(names(Hs578T_df_ls), function(x){
  df <- Hs578T_df_ls[[x]]
  out <- ggplot(df, aes(x = log2FC, y = mlog10FDR, color = signif)) + geom_point_rast(size = 0.8) + ArchR::theme_ArchR() +
    geom_hline(yintercept = -log10(0.1), lty = "dashed") + geom_vline(xintercept = c(0), lty = "dashed") + 
    scale_color_manual(values = c("N" = "lightgray", "siFOXD1_Up" = "red", "siFOXD1_Dn" = "blue")) +
    ggtitle(x) + labs(x = "log2FC", y = "-log10(FDR)")
  return(out)
})

pdf("output/Plots/v2_C_Hs578T_VolcanoPlot.pdf", width = 5, height = 5)
p6
dev.off()

#export
peaks <- rowRanges(Hs578T_se)
gr2 <- GRangesList(Hs578T_siFOXD1_1_Up = peaks[rownames(Hs578T_df_ls$Hs578T_siFOXD1_1)[Hs578T_df_ls$Hs578T_siFOXD1_1$signif == "siFOXD1_Up"]], 
                   Hs578T_siFOXD1_1_Dn = peaks[rownames(Hs578T_df_ls$Hs578T_siFOXD1_1)[Hs578T_df_ls$Hs578T_siFOXD1_1$signif == "siFOXD1_Dn"]], 
                   Hs578T_siFOXD1_2_Up = peaks[rownames(Hs578T_df_ls$Hs578T_siFOXD1_2)[Hs578T_df_ls$Hs578T_siFOXD1_2$signif == "siFOXD1_Up"]], 
                   Hs578T_siFOXD1_2_Dn = peaks[rownames(Hs578T_df_ls$Hs578T_siFOXD1_2)[Hs578T_df_ls$Hs578T_siFOXD1_2$signif == "siFOXD1_Dn"]])
lapply(names(gr2), function(x){
  g <- gr2[[x]]
  d <- data.frame(seqnames = seqnames(g), start = start(g)-1, end = end(g))
  write.table(d, paste0("output/output_bed_v2/", x, ".bed"), row.names = F, col.names = F, quote = F, sep = "\t")
})

#----------------- chromVAR -----------------#
library(chromVAR)
library(chromVARmotifs)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

#----- BT549
#add GC bias
BT549_se <- addGCBias(BT549_se, genome = BSgenome.Hsapiens.UCSC.hg38)
#homer motif
motifs <- chromVARmotifs::human_pwms_v2
motif_ix <- matchMotifs(motifs, BT549_se, genome = BSgenome.Hsapiens.UCSC.hg38)
#calculate deviations
bg <- getBackgroundPeaks(object = BT549_se)
dev <- computeDeviations(object = BT549_se, annotations = motif_ix, background_peaks = bg)

BT549_deltaMS <- rowMeans(assays(dev)$z[,c(5:6)]) - rowMeans(assays(dev)$z[,c(1:4)])
names(BT549_deltaMS) <- rowData(dev)$name
sort(BT549_deltaMS, decreasing = T)

#----- Hs578T
#add GC bias
Hs578T_se <- addGCBias(Hs578T_se, genome = BSgenome.Hsapiens.UCSC.hg38)
#homer motif
motif_ix2 <- matchMotifs(motifs, Hs578T_se, genome = BSgenome.Hsapiens.UCSC.hg38)
#calculate deviations
bg2 <- getBackgroundPeaks(object = Hs578T_se)
dev2 <- computeDeviations(object = Hs578T_se, annotations = motif_ix2, background_peaks = bg2)

Hs578T_deltaMS <- rowMeans(assays(dev2)$z[,c(5:6)]) - rowMeans(assays(dev2)$z[,c(1:4)])
names(Hs578T_deltaMS) <- rowData(dev2)$name
sort(Hs578T_deltaMS, decreasing = T)

#----- merge
idx1 <- head(sort(BT549_deltaMS, decreasing = T),30) %>% names
idx2 <- head(sort(Hs578T_deltaMS, decreasing = T),30) %>% names
top30 <- intersect(idx1, idx2)

df <- data.frame(BT549 = BT549_deltaMS, Hs578T = Hs578T_deltaMS, label = names(BT549_deltaMS), col = "gray")

df$col[which(df$label %in% idx1)] <- "orange"
df$col[which(df$label %in% idx2)] <- "purple"
df$col[which(df$label %in% top10)] <- "red"
df$label[which(df$label %ni% c(top10, "FOXD1"))] <- ""

p7 <- ggplot(df, aes(x = BT549, y = Hs578T, fill = col, label = label)) + geom_point(size = 2, pch = 21) + ArchR::theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + geom_hline(yintercept = 0, lty = "dashed") +
  labs(x = "delta ChromVAR scores, siNC - siFOXD1 in BT549", y = "delta ChromVAR scores, siNC - siFOXD1 in Hs578T") +
  ggrepel::geom_label_repel(max.overlaps = 100) + scale_fill_identity()

pdf("output/Plots/v2_C_ChromVARscore_Scatter.pdf", width = 11, height = 11)
p7
dev.off()
pdf("output/Plots/v2_C_ChromVARscore_Scatter_2.pdf", width = 5, height = 5)
p7
dev.off()

write.csv(df, "output/Tables/v2_ChromVARmotifScore.csv")

#----- GO enrichment visualization -----#
go_list <- list.files("output/great_output/", pattern = "GREAT_")
GREAT_go <- lapply(go_list, function(i){
  out <- data.table::fread(paste0("output/great_output/", i), header = F, skip = 1)[c(1:20), c(1,4)] %>% data.frame %>% `colnames<-`(., c("Term", "FDR"))
  out$mlog10FDR <- -log10(out$FDR)
  return(out)
})
names(GREAT_go) <- go_list

p9 <- lapply(names(GREAT_go), function(i){
  out <- ggplot(GREAT_go[[i]], aes(x = mlog10FDR, y = reorder(Term, mlog10FDR))) +
    geom_bar(stat = "identity") + ArchR::theme_ArchR() + ggtitle(i)
  return(out)})

pdf("output/Plots/v2_C_GREAT.pdf", height = 5, width = 10)
p9
dev.off()
