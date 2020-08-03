## Run dmrseq
# oct 21 2019

suppressPackageStartupMessages({
  library(dmrseq)
  library(dplyr)
  library(ggplot2)
  library(annotatr)
  library(BiocParallel)
})

source("scripts/get_table_with_annots.R")
metafile <- "data/metadata_crc.txt"
#print(metafile)
folder <- "/home/Shared_s3it/sorjuela/CRC_outputs"


#Read in metadata
metadata <- read.delim(metafile, header = TRUE, as.is = TRUE, sep = "\t")
infile <- file.path(folder, "methextract",sprintf("%s_pe.CpG_report.txt", metadata$names))
#infile2 <- sprintf("data/methextract/%s_pe.bismark.cov.gz", metadata$names)

#metadata$patient <- gsub("[0-9]+\\.A-", "", metadata$names)
#metadata$patient <- gsub("_[a-z]+","",metadata$patient)

metadata$age_group <- ifelse(metadata$age > 40, "old", "young")


#Read in cov files
bismarkBSseq <- read.bismark(files = infile,
                             rmZeroCov = TRUE,
                             strandCollapse = TRUE,
                             verbose = TRUE, 
                             colData = metadata,
                             BPPARAM = MulticoreParam(1),
                             nThread = 1,
                             BACKEND = "HDF5Array",
                             dir = "data/HDF5_cov_crc",
                             replace = TRUE) 

save(bismarkBSseq, file = "data/bsseq_cpgreport.RData")
#load("data/bsseq_cpgreport.RData")

#Filter coverage low
loci.idx <- DelayedMatrixStats::rowSums2(getCoverage(bismarkBSseq, type="Cov") >= 10 ) >= 10
bismarkBSseq <- bismarkBSseq[loci.idx,] 

#Filter chrom X,Y and weird ones
loci.idx <- seqnames(bismarkBSseq) %in% c(1:22)
bismarkBSseq <- bismarkBSseq[loci.idx,] #2,853,792

#Filter/reduce coverage high (test)
#loci.idx <- DelayedMatrixStats::rowSums2(getCoverage(bismarkBSseq, type="Cov") <= 200 ) == 36
#bismarkBSseq <- bismarkBSseq[loci.idx,] 

#or from BiSeq
cov <- getCoverage(bismarkBSseq, type = "Cov")
meth <- getCoverage(bismarkBSseq, type = "M")
ind.cov <- cov > 0
quant <- quantile(cov[ind.cov], 0.9)
quant  

limitCov <- function(cov, meth, maxCov, object){
  cov <- as.matrix(cov)
  meth <- as.matrix(meth)
  indCov <- cov > maxCov
  fraction <- meth[indCov] / cov[indCov]
  cov[indCov] <- as.integer(maxCov)
  meth[indCov] <- as.integer(round(fraction * maxCov))
  object <- BSseq(M = meth, 
                  Cov = cov, 
                  pData = colData(object), 
                  gr = granges(object),
                  sampleNames = colData(bismarkBSseq)$names)
  return(object)
}

bismarkBSseq <- limitCov(cov, meth, quant, bismarkBSseq)
colData(bismarkBSseq)$age_group <- as.factor(colData(bismarkBSseq)$age_group)

#Generate bws
cov <- getCoverage(bismarkBSseq, type = "Cov")
meth <- getCoverage(bismarkBSseq, type = "M")
meth_vals <- meth /cov * 100
colnames(meth_vals) <- colData(bismarkBSseq)$names

make_bigwigs <- function(sample, bs, meth_vals, folder, chromsizesFile){
  #make dataframe with midpoint and score for file
  bg <- data.frame(chr = GenomeInfoDb::seqnames(bs),
                   start = start(bs),
                   end = end(bs),
                   score = abs(meth_vals[,grep(sample, colnames(meth_vals))]),
                   stringsAsFactors = F,
                   row.names = NULL)
  
  #remove NAs (this is just to double check, the user should input a filtered
  #matrix)
  bg <- bg[!is.na(meth_vals[,grep(sample, colnames(meth_vals))]),]
  
  if(folder == "."){
    folder = ""
  } else {folder = paste0(folder,"/")}
  
  #make bedgraph
  scorename <- "meth"
  write.table(bg, file = sprintf("%s%s.%s.bedgraph", folder, scorename, sample),
              sep = "\t", row.names = F, col.names = F, quote = F)
  
  #sort file
  cmd1 <- sprintf(
    "sort -k1,1 -k2,2n %s%s.%s.bedgraph > %s%s.%s.sorted.bedgraph",
    folder, scorename, sample, folder, scorename, sample)
  cat(cmd1, "\n")
  system(cmd1)
  
  #then make bigwig
  cmd2 <- sprintf("bedGraphToBigWig %s%s.%s.sorted.bedgraph %s %s%s.%s.bw",
                  folder, scorename, sample, chromsizesFile, folder, 
                  scorename, sample)
  cat(cmd2, "\n")
  system(cmd2)
}

sapply(colnames(meth_vals), make_bigwigs, 
       bs = bismarkBSseq, 
       meth_vals = meth_vals,
       folder = ".", 
       chromsizesFile = "/home/Shared_taupo/stephany/reference/hg19.chrom.sizes.modGL")


#Get diagnostics
mean_meth <- rowMeans(meth_vals) 
mean_cov <- rowMeans(cov)
var_meth <- rowVars(as.matrix(meth_vals))

d <- data.frame(mean_meth = mean_meth,
                mean_cov =  mean_cov,
                #diff = diffs,
                var_meth = var_meth)

#mean meth Vs mean cov
ggplot(d) + aes(x=mean_cov, y=mean_meth) + 
  geom_bin2d() + 
  scale_fill_distiller(palette='RdBu', trans='log10') + 
  scale_x_continuous(breaks=scales::pretty_breaks()) +
  theme_bw()
ggsave("figures/covVsmeth_crc.png")

#MDS
meth_vals <- meth /cov
methsTR <- asin(2*meth_vals-1)
#bad <- rowSums(is.finite(methsTR)) < ncol(methsTR)
#if(any(bad)) methsTR <- methsTR[!bad,,drop=FALSE]

mds_meth <- limma::plotMDS(methsTR, top = 10000, plot = FALSE)$cmdscale.out

df <- data.frame(dim1 = mds_meth[,1], dim2 = mds_meth[,2],
                 names = colnames(meth_vals), treat = metadata$age_group)

ggplot()+
  geom_point(data = df, mapping = aes_(x=~dim1, y=~dim2, color=~treat)) +
  geom_text(data = df, mapping = aes_(x=~dim1, y=~dim2-0.01, label=~names)) +
  theme_bw()
ggsave("figures/MDS_top10mil_crc.png")

#MV
ggplot(d) + aes(x=mean_cov, y=var_meth) + 
  geom_bin2d() + 
  scale_fill_distiller(palette='RdBu', trans='log10') + 
  scale_x_continuous(breaks=scales::pretty_breaks()) +
  theme_bw()
ggsave("figures/meanVsvar_crc.png")

#meth per sample
pos <- paste0(seqnames(bismarkBSseq),".", start(bismarkBSseq))
meth_vals <- as.data.frame(meth_vals)
meth_vals$pos <- pos
meth_vals_melt <- reshape2::melt(meth_vals, id.vars = "pos", measure.vars = metadata$names)

ggplot(meth_vals_melt) + geom_violin(aes(x = variable, y = value)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("figures/meth_violins_crc.png")


#### Run dmrseq ####
#scale beta func

#nonCIMP

bsn <- bismarkBSseq[,grepl("non",metadata$lesion)]
set.seed(1234)
DMRsage_non <- dmrseq(bs=bsn,
                        testCovariate="lesion", 
                        cutoff = 0.05, 
                        BPPARAM = MulticoreParam(1),
                        adjustCovariate = "patient",
                        maxGap = 100,
                        maxGapSmooth = 1000,
                        matchCovariate = "gender",
                        minNumRegion = 3)


DMRs_annot <- get_table_with_annots(DMRsage_non)
DMRs_annot$beta <- -DMRs_annot$beta
save(DMRs_annot,file = "data/DMRs_nonCIMP.RData")

df <- as.data.frame(DMRs_annot)
write.table(df, file = "data/DMRs_nonCIMP.txt", quote = FALSE, row.names = FALSE)

#CIMP
bsc <- bismarkBSseq[,!grepl("non",metadata$lesion)]
set.seed(1234)
DMRs_cimp <- dmrseq(bs=bsc,
                      testCovariate="lesion", 
                      cutoff = 0.05, 
                      BPPARAM = MulticoreParam(1),
                      adjustCovariate = "patient",
                      maxGap = 100,
                      maxGapSmooth = 1000,
                      matchCovariate = "gender",
                      minNumRegion = 3)


DMRs_annot <- get_table_with_annots(DMRs_cimp)
DMRs_annot$beta <- -DMRs_annot$beta
save(DMRs_annot,file = "data/DMRs_CIMP.RData")

df <- as.data.frame(DMRs_annot)
write.table(df, file = "data/DMRs_CIMP.txt", quote = FALSE, row.names = FALSE)

#both CRCs
colData(bismarkBSseq)$tissue <- ifelse(grepl("Normal", bismarkBSseq$lesion), 
                                       "Normal", "CRC")
set.seed(1234)
DMRs_cimp <- dmrseq(bs=bismarkBSseq,
                    testCovariate="tissue", 
                    cutoff = 0.05, 
                    BPPARAM = MulticoreParam(4),
                    adjustCovariate = "patient",
                    maxGap = 100,
                    maxGapSmooth = 1000,
                    matchCovariate = "gender",
                    minNumRegion = 3)


DMRs_annot <- get_table_with_annots(DMRs_cimp)
DMRs_annot$beta <- -DMRs_annot$beta
save(DMRs_annot,file = "data/DMRs_allCRC.RData")

df <- as.data.frame(DMRs_annot)
write.table(df, file = "data/DMRs_allCRC.txt", quote = FALSE, row.names = FALSE)