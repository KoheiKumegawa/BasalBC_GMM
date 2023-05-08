#----------------------------------------------------------------------------
# TCGAdata_analysis.R
#----------------------------------------------------------------------------
library(SummarizedExperiment)
library(mclust)
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(ggridges)
library(ComplexHeatmap)
library(circlize)
se <- readRDS("rds/TCGA_BRCA_RNAExp.rds")
assays(se)$fpkm_unstrand_log2 <- log2(assays(se)$fpkm_unstrand+1)

#------------------- basal -------------------#
#primary, basal tumors, 194 samples
basal_se <- se[, which(se$sample_type == "Primary Tumor" & se$paper_BRCA_Subtype_PAM50 == "Basal")]
expStats <- data.frame(mean = rowMeans(assays(basal_se)$fpkm_unstrand_log2),
                       var  = rowVars(assays(basal_se)$fpkm_unstrand_log2),
                       symbol = mcols(basal_se)$gene_name)
p1 <- ggplot(expStats, aes(x = var, y = mean)) + geom_hex(bins = 100) + 
  ArchR::theme_ArchR() + scale_fill_viridis_c(trans = "log", option = "E") + 
  labs(x = "Variance", y = "Mean Expression") + 
  geom_hline(yintercept = 0.5, linetype = "dotted") + 
  geom_vline(xintercept = 0.5, linetype = "dotted")
selected_genes <- rownames(expStats)[which(expStats$mean > 0.5 & expStats$var > 0.5)]
basal_se <- basal_se[selected_genes,]

#mclust
mclust1 <- parallel::mclapply(rownames(basal_se), function(x){
  out <- Mclust(assays(basal_se)$fpkm_unstrand_log2[x,], modelNames = "V") 
  return(out)
}, mc.cores = 16)
names(mclust1) <- mcols(basal_se)$gene_name

idx <- lapply(seq_along(mclust1), function(x) length(unique(mclust1[[x]]$classification))) %>% unlist()
p2 <- ggplot(data.frame(table(idx)), aes(x = "", y = Freq, fill = idx)) + geom_bar(stat="identity", width = 1, alpha = 0.8) + 
      coord_polar("y", start=0) + ArchR::theme_ArchR()

mclust2 <- mclust1[which(idx != 1)]

#at least 10 samples in each groups
idy <- lapply(seq_along(mclust2), function(x){
  a <- table(mclust2[[x]]$classification)
  out <- min(a) >= 10
  return(out)
}) %>% unlist()
mclust3 <- mclust2[idy]

#survival analysis
survP <- parallel::mclapply(seq_along(mclust3), function(x){
  group <- mclust3[[x]]$classification
  group <- group[colnames(basal_se)]
  surv <- data.frame(sample = colnames(basal_se),
                     daysFO = basal_se$days_to_last_follow_up,
                     daysDE = basal_se$days_to_death,
                     status = basal_se$vital_status,
                     group  = group)
  surv$days <- surv$daysFO
  surv$days[is.na(surv$days)] <- surv$daysDE[is.na(surv$days)]
  surv$status2 <- 0
  surv$status2[surv$status == "Dead"] <- 1 
  
  sd <- survdiff(Surv(surv$days, surv$status2)~surv$group)
  p <- pchisq(sd$chisq, length(sd$n)-1, lower.tail = FALSE)
  return(p)
}, mc.cores = 16)
survP <- unlist(survP)

survPDF <- data.frame(name = names(mclust3),
                      LogRankP = survP,
                      LogRankPlog = -log10(survP))
survPDF <- survPDF[order(survPDF$LogRankP),]
survPDF$rank <- c(1:nrow(survPDF))

write.csv(survPDF, "output/Tables/S1_basal_logRankTest.csv")
p3 <- ggplot(survPDF, aes(x = rank, y = LogRankPlog)) + geom_point() + ArchR::theme_ArchR() +
      geom_hline(yintercept = -log10(0.01), linetype = "dashed") + labs(x = "Rank-sorted genes", y = "-log2(P-value)")

#ridge plot
sig_genes <- survPDF$name[survPDF$LogRankPlog >= 2]
p4 <- lapply(sig_genes, function(x){
  idx <- rownames(basal_se)[mcols(basal_se)$gene_name == x]
  group <- mclust3[[x]]$classification
  df <- data.frame(sample = colnames(basal_se),
                   group  = paste0("group", group[colnames(basal_se)]),
                   exp    = assays(basal_se)$fpkm_unstrand_log2[idx,])
  p <- ggplot(df, aes(x = exp, y = group, fill = group)) + geom_density_ridges(alpha = 0.5) + 
       ArchR::theme_ArchR() + ggtitle(x)
  return(p)
})

#survival plot
p5 <- list()
for(i in seq_along(sig_genes)){
  x <- sig_genes[i]
  group <- mclust3[[x]]$classification
  surv <- data.frame(sample = colnames(basal_se),
                     daysFO = basal_se$days_to_last_follow_up,
                     daysDE = basal_se$days_to_death,
                     status = basal_se$vital_status,
                     group  = group[colnames(basal_se)])
  surv$days <- surv$daysFO
  surv$days[is.na(surv$days)] <- surv$daysDE[is.na(surv$days)]
  surv$status2 <- 0
  surv$status2[surv$status == "Dead"] <- 1 
  sf <- survfit(Surv(surv$days, surv$status2)~surv$group)
  p <- ggsurvplot(fit = sf, data = surv,
                  pval = TRUE, pval.method = TRUE,
                  risk.table = TRUE, conf.int = FALSE,
                  ncensor.plot = FALSE, size = 1.5,
                  title = x,
                  legend.title = "Expression Pattern",
                  log.rank.weights = "1",
                  risk.table.title = "",
                  risk.table.y.text.col = TRUE,
                  risk.table.y.text = FALSE)
  p5[[i]] <- p
}

#output
pdf("output/Plots/S1_RNA_Variance.pdf", width = 3.5, height = 4)
p1
dev.off()
pdf("output/Plots/S1_No_Mclusters.pdf", width = 4, height = 4)
p2
dev.off()
pdf("output/Plots/S1_LogRankTest.pdf", width = 4, height = 4)
p3
dev.off()
pdf("output/Plots/S1_ExpRidges.pdf", width = 4, height = 4)
p4
dev.off()
pdf("output/Plots/S1_Survplot.pdf", width = 6, height = 8)
p5
dev.off()

#------------------- FOXD1 in other subtypes -------------------#
other_se <- list(LumA = se[, which(se$sample_type == "Primary Tumor" & se$paper_BRCA_Subtype_PAM50 == "LumA")],
                 LumB = se[, which(se$sample_type == "Primary Tumor" & se$paper_BRCA_Subtype_PAM50 == "LumB")],
                 Her2 = se[, which(se$sample_type == "Primary Tumor" & se$paper_BRCA_Subtype_PAM50 == "Her2")],
                 Normal = se[, which(se$sample_type == "Primary Tumor" & se$paper_BRCA_Subtype_PAM50 == "Normal")])
mclustFOXD1 <- parallel::mclapply(other_se, function(x){
  out <- Mclust(assays(x)$fpkm_unstrand_log2["ENSG00000251493.5",], modelNames = "V") 
  return(out)
}, mc.cores = 16)

p6 <- parallel::mclapply(names(other_se), function(x){
  se <- other_se[[x]]
  group <- mclustFOXD1[[x]]$classification
  df <- data.frame(sample = colnames(se),
                   group  = paste0("group", group[colnames(se)]),
                   exp    = assays(se)$fpkm_unstrand_log2["ENSG00000251493.5",])
  p <- ggplot(df, aes(x = exp, y = group, fill = group)) + geom_density_ridges(alpha = 0.5) + ArchR::theme_ArchR() + ggtitle(x)
  return(p)
})

p7 <- list()
for(i in seq_along(other_se)){
  x <- names(other_se)[i]
  tmp_se <- other_se[[x]]
  group <- mclustFOXD1[[x]]$classification
  surv <- data.frame(sample = colnames(tmp_se),
                     daysFO = tmp_se$days_to_last_follow_up,
                     daysDE = tmp_se$days_to_death,
                     status = tmp_se$vital_status,
                     group  = group[colnames(tmp_se)])
  surv$days <- surv$daysFO
  surv$days[is.na(surv$days)] <- surv$daysDE[is.na(surv$days)]
  surv$status2 <- 0
  surv$status2[surv$status == "Dead"] <- 1 
  sf <- survfit(Surv(surv$days, surv$status2)~surv$group)
  p <- ggsurvplot(fit = sf, data = surv,
                  pval = TRUE, pval.method = TRUE,
                  risk.table = TRUE, conf.int = FALSE,
                  ncensor.plot = FALSE, size = 1.5, #linetype = c(1, 3),
                  title = x,
                  legend.title = "Expression Pattern",
                  log.rank.weights = "1",
                  risk.table.title = "",
                  risk.table.y.text.col = TRUE,
                  risk.table.y.text = FALSE)
  p7[[i]] <- p
}

pdf("output/Plots/S2_ExpRidge_FOXD1_otherSubtypes.pdf", width = 4, height = 4)
p6
dev.off()
pdf("output/Plots/S2_Survplot_FOXD1_otherSubtypes.pdf", width = 6, height = 8)
p7
dev.off()
