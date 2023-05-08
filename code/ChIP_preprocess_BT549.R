#------------------------------------------------------------------------------
# ChIP_preprocess_BT549.R
#------------------------------------------------------------------------------
library(Rsamtools)
library(data.table)
library(dplyr)
library(rtracklayer)
library(SummarizedExperiment)
'%ni%' <- Negate("%in%")

#--------- functions ---------#
bamToFragmentGR <- function(
  bamPATH = NULL,
  bamNAME = NULL,
  offsetPlus = 4,
  offsetMinus = -5,
  bamFlag = NULL
){
  if(is.null(bamPATH)){
    stop("Please set PATH to bam files")
  }
  if(is.null(bamNAME)){
    stop("No input bamNAME; please recheck your input")
  }
  if(is.null(bamFlag)){
    stop("Please set bamFlag using Rsamtools's scanBamFlag!")
  }
  
  #1. Read In Bam File
  sF <- scanBam(bamPATH, param = ScanBamParam(flag = bamFlag, what = c("rname","pos", "isize")))[[1]]
  
  #2. Make Fragment Table
  dt <- data.table(seqnames = sF$rname, start = sF$pos + offsetPlus, end = sF$pos + abs(sF$isize) - 1 + offsetMinus)
  
  #3. Make fragment Granges and remove unwanted chromosomes
  gr <- GRanges(seqnames = dt$seqnames, IRanges(start = dt$start, end = dt$end))
  idy = which(seqnames(gr) %in% seqlevels(gr)[grep("random|chrM|chrUn|chrEBV", seqlevels(gr))])
  gr <- gr[-idy]
  gr <- dropSeqlevels(gr, seqlevels(gr)[grep("random|chrM|chrUn|chrEBV", seqlevels(gr))])
  mcols(gr) <- DataFrame(sample = bamNAME)
  
  #4. output Granges List
  return(gr)
}

MakeSamplePeakSet <- function(gr, by = "score"){
  #nonOverlappingGRanges
  stopifnot(by %in% colnames(mcols(gr)))
  
  #function for picking up most significant peaks
  clusterGRanges <- function(gr, by = "score"){
    gr <- sort(sortSeqlevels(gr))
    r <- GenomicRanges::reduce(gr, min.gapwidth=0L, ignore.strand=TRUE)
    o <- findOverlaps(gr,r)
    mcols(gr)$cluster <- subjectHits(o)
    gr <- gr[order(mcols(gr)[,by], decreasing = TRUE),]
    gr <- gr[!duplicated(mcols(gr)$cluster),]
    gr <- sort(sortSeqlevels(gr))
    mcols(gr)$cluster <- NULL
    return(gr)
  }
  
  #iteration of filtering overlapping peaks
  i <-  0
  gr_converge <- gr
  while(length(gr_converge) > 0){
    i <-  i + 1
    gr_selected <- clusterGRanges(gr = gr_converge, by = by)
    gr_converge <- subsetByOverlaps(gr_converge, gr_selected, invert=TRUE) #blacklist selected gr
    if(i == 1){ #if i=1 then set gr_all to clustered
      gr_all <- gr_selected
    }else{
      gr_all <- c(gr_all, gr_selected)
    }
  }
  gr_all <- sort(sortSeqlevels(gr_all))
  return(gr_all)
}

#--------- analysis ---------#
bamList <- list.files("data/", pattern = "BT549") %>% grep("bam$", ., value = T)
names(bamList) <- c("BT549-siFOXD1-1-1","BT549-siFOXD1-1-2","BT549-siFOXD1-2-1","BT549-siFOXD1-2-2","BT549-siNC-1-1","BT549-siNC-1-2")

#bam to gr
fragments <- parallel::mclapply(names(bamList), 
                                function(x){bamToFragmentGR(bamPATH  = paste0("data/", bamList[[x]]),
                                                            bamNAME = x,
                                                            bamFlag = scanBamFlag(isMinusStrand = FALSE, isProperPair  = TRUE)
                                )}, mc.cores = 6)
names(fragments) <- names(bamList)

#Make sample peak set
summitPeaks <- list.files("data/", pattern = c("BT549")) %>% grep("summits", ., value = T)
gr_ls <- lapply(summitPeaks, function(x) import(paste0("data/", x)))
names(gr_ls) <- names(bamList)

BSgenome   <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
chromSizes <- GRanges(names(seqlengths(BSgenome)), IRanges(1, seqlengths(BSgenome))) %>% 
  GenomeInfoDb::keepStandardChromosomes(., pruning.mode = "coarse")
blacklist  <- rtracklayer::import.bed("ref/hg38-blacklist.v2.bed")

gr_ls_proc <- parallel::mclapply(gr_ls, function(x){
  gr <- resize(x, width = 501, fix = "center") %>%
    subsetByOverlaps(., chromSizes, type = "within") %>%
    subsetByOverlaps(., blacklist, invert=TRUE) %>%
    MakeSamplePeakSet(., by = "score")
  mcols(gr)$scorePerMillion <- mcols(gr)$score / (sum(mcols(gr)$score) / 1000000)
  return(gr)
}, mc.cores = 10) %>% GenomicRanges::GRangesList(.)

#make consensus peakset
gr_cumulative <- MakeSamplePeakSet(unlist(gr_ls_proc), by = "scorePerMillion")
mcols(gr_cumulative)$sampleOverlap <- countOverlaps(gr_cumulative, gr_ls_proc)
reproduciblePeaks <- gr_cumulative[which(mcols(gr_cumulative)$scorePerMillion >= 5 &
                                           mcols(gr_cumulative)$sampleOverlap >= 2 &
                                           seqnames(gr_cumulative) %ni% c("chrY", "chrM"))]
names(reproduciblePeaks) <- paste0("BT549_K27AC_", c(1:length(reproduciblePeaks)))
mcols(reproduciblePeaks)$name <- names(reproduciblePeaks)


#make SE
fragments <- unlist(as(fragments, "GRangesList"))
overlapDF <- DataFrame(findOverlaps(reproduciblePeaks, fragments, ignore.strand = TRUE, maxgap=-1L, minoverlap=0L, type = "any"))
overlapDF$name <- mcols(fragments)[overlapDF[, 2], "sample"]
overlapTDF <- transform(overlapDF, id = match(name, unique(name)))
#Summarize
sparseM <- Matrix::sparseMatrix(
  i = overlapTDF[, 1], 
  j = overlapTDF[, 4],
  x = rep(1, nrow(overlapTDF)), 
  dims = c(length(reproduciblePeaks), length(unique(overlapDF$name))))
colnames(sparseM) <- unique(overlapDF$name)
rownames(sparseM) <- reproduciblePeaks$name
sparseM.cpm <- edgeR::cpm(sparseM, log = TRUE, prior.count = 1)
rownames(sparseM.cpm) <- reproduciblePeaks$name
#SummarizedExperiment
se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = sparseM, log2cpm = sparseM.cpm),rowRanges = reproduciblePeaks)

saveRDS(se, "rds/BT549_K27AC_se_v2.rds")
