#------------------------------------------------------------------------------
# RNA_analysis.R
#------------------------------------------------------------------------------
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(ggrastr)
source("code/edgeR_PairwiseFunction.R")

se <- readRDS("rds/se.rds")

#differential analysis
diff_ls <- list(BT549_siFOXD1_1 = list("BT549_siFOXD1_1", "BT549_siNC_1"), 
                BT549_siFOXD1_2 = list("BT549_siFOXD1_2","BT549_siNC_1"), 
                Hs578T_siFOXD1_1 = list("Hs578T_siFOXD1_1","Hs578T_siNC_1"),
                Hs578T_siFOXD1_2 = list("Hs578T_siFOXD1_2","Hs578T_siNC_1"))
DiffTest <- lapply(diff_ls, function(x) edgeR_pairwise(se, compareCol = "sampleType", topGroup = x[[1]], bottomGroup = x[[2]]))
DiffTest


df_ls <- lapply(names(DiffTest), function(x){
  df <- assay(DiffTest[[x]])[, c("log2FC", "FDR")] %>% data.frame
  df$mlog10FDR <- -log10(df$FDR)
  df$signif <- "N"
  df$signif[which(df$log2FC > 1 & df$FDR < 0.01)] <- "siFOXD1_Up"
  df$signif[which(df$log2FC < -1 & df$FDR < 0.01)] <- "siFOXD1_Dn"
  df$signif <- factor(df$signif, levels = c("N", "siFOXD1_Up", "siFOXD1_Dn"))
  df <- df[order(df$signif),]
  return(df)
})
names(df_ls) <- names(DiffTest)
                  
p1 <- lapply(names(df_ls), function(x){
  df <- df_ls[[x]]
  out <- ggplot(df, aes(x = log2FC, y = mlog10FDR, color = signif)) + geom_point_rast(size = 0.8) + ArchR::theme_ArchR() +
    geom_hline(yintercept = 2, lty = "dashed") + geom_vline(xintercept = c(-1, 1), lty = "dashed") + 
    scale_color_manual(values = c("N" = "lightgray", "siFOXD1_Up" = "red", "siFOXD1_Dn" = "blue")) +
    ggtitle(x) + labs(x = "log2FC", y = "-log10(FDR)") + ylim(0, max(df$mlog10FDR)*0.9) + xlim(-8,8)
  return(out)
})

pdf("output/Plots/RNA_VolcanoPlot.pdf", width = 5, height = 5)
p1
dev.off()

lapply(names(df_ls), function(i) write.csv(df_ls[[i]], paste0("output/Tables/DEG_", i, ".csv")))

BT549_siFOXD1_1_Up <- rownames(df_ls$BT549_siFOXD1_1)[which(df_ls$BT549_siFOXD1_1$signif == "siFOXD1_Up")]
BT549_siFOXD1_2_Up <- rownames(df_ls$BT549_siFOXD1_2)[which(df_ls$BT549_siFOXD1_2$signif == "siFOXD1_Up")]
BT549_siFOXD1_Up <- intersect(BT549_siFOXD1_1_Up, BT549_siFOXD1_2_Up)

BT549_siFOXD1_1_Dn <- rownames(df_ls$BT549_siFOXD1_1)[which(df_ls$BT549_siFOXD1_1$signif == "siFOXD1_Dn")]
BT549_siFOXD1_2_Dn <- rownames(df_ls$BT549_siFOXD1_2)[which(df_ls$BT549_siFOXD1_2$signif == "siFOXD1_Dn")]
BT549_siFOXD1_Dn <- intersect(BT549_siFOXD1_1_Dn, BT549_siFOXD1_2_Dn)

Hs578T_siFOXD1_1_Up <- rownames(df_ls$Hs578T_siFOXD1_1)[which(df_ls$Hs578T_siFOXD1_1$signif == "siFOXD1_Up")]
Hs578T_siFOXD1_2_Up <- rownames(df_ls$Hs578T_siFOXD1_2)[which(df_ls$Hs578T_siFOXD1_2$signif == "siFOXD1_Up")]
Hs578T_siFOXD1_Up <- intersect(Hs578T_siFOXD1_1_Up, Hs578T_siFOXD1_2_Up)

Hs578T_siFOXD1_1_Dn <- rownames(df_ls$Hs578T_siFOXD1_1)[which(df_ls$Hs578T_siFOXD1_1$signif == "siFOXD1_Dn")]
Hs578T_siFOXD1_2_Dn <- rownames(df_ls$Hs578T_siFOXD1_2)[which(df_ls$Hs578T_siFOXD1_2$signif == "siFOXD1_Dn")]
Hs578T_siFOXD1_Dn <- intersect(Hs578T_siFOXD1_1_Dn, Hs578T_siFOXD1_2_Dn)

write.table(BT549_siFOXD1_Up, "output/Tables/DEG_BT549_siFOXD1_Up.txt", quote = F, row.names = F, col.names = F)
write.table(BT549_siFOXD1_Dn, "output/Tables/DEG_BT549_siFOXD1_Dn.txt", quote = F, row.names = F, col.names = F)
write.table(Hs578T_siFOXD1_Up, "output/Tables/DEG_Hs578T_siFOXD1_Up.txt", quote = F, row.names = F, col.names = F)
write.table(Hs578T_siFOXD1_Dn, "output/Tables/DEG_Hs578T_siFOXD1_Dn.txt", quote = F, row.names = F, col.names = F)
