################################################################################
# R script to read-in CpG-report files from bismark (in snakemake workflow), 
# limit coverage, generate diagnostic plots, and detect differentially methylated
# regions in healthy samples
#
# oct 21 2019
################################################################################

suppressPackageStartupMessages({
  library(dmrseq)
  library(dplyr)
  library(ggplot2)
  library(annotatr)
  library(BiocParallel)
})

source("scripts/get_table_with_annots.R")
metafile <- "metadata_healthy.txt"
#print(metafile)

#Read in metadata
metadata <- read.delim(metafile, header = TRUE, as.is = TRUE, sep = "\t")
infile <- sprintf("data/methextract/%s_pe.CpG_report.txt", metadata$names)

metadata$age_group <- ifelse(metadata$age > 50, "old", "young")


#Read in cov files
bismarkBSseq <- read.bismark(files = infile,
                             rmZeroCov = TRUE,
                             strandCollapse = TRUE,
                             verbose = TRUE, 
                             colData = metadata,
                             BPPARAM = MulticoreParam(1),
                             nThread = 1,
                             BACKEND = "HDF5Array",
                             dir = "data/HDF5_cov",
                             replace = TRUE) #26,660,520

#save(bismarkBSseq, file = "data/rdata/bsseq_cpgreport.RData")

#Filter coverage low
loci.idx <- DelayedMatrixStats::rowSums2(getCoverage(bismarkBSseq, type="Cov") >= 10 ) >= 10
bismarkBSseq <- bismarkBSseq[loci.idx,] #2,481,876

#Filter chrom X,Y and weird ones
loci.idx <- seqnames(bismarkBSseq) %in% c(1:22)
bismarkBSseq <- bismarkBSseq[loci.idx,] #2,418,647

#limit coverage
cov <- getCoverage(bismarkBSseq, type = "Cov")
meth <- getCoverage(bismarkBSseq, type = "M")
ind.cov <- cov > 0
quant <- quantile(cov[ind.cov], 0.9)
quant  #66

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
#save(bismarkBSseq, file="data/rdata/bsseq_age_filt.RData") #too large for github

#### Diagnostic plots ####
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
ggsave("figures/covVsmeth.png")

#MDS
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
ggsave("figures/MDS_top10mil.png")

#MV
ggplot(d) + aes(x=mean_cov, y=var_meth) + 
  geom_bin2d() + 
  scale_fill_distiller(palette='RdBu', trans='log10') + 
  scale_x_continuous(breaks=scales::pretty_breaks()) +
  theme_bw()
ggsave("figures/meanVsvar.png")

#meth per sample
pos <- paste0(seqnames(bismarkBSseq),".", start(bismarkBSseq))
meth_vals <- as.data.frame(meth_vals)
meth_vals$pos <- pos
meth_vals_melt <- reshape2::melt(meth_vals, id.vars = "pos", measure.vars = metadata$names)

ggplot(meth_vals_melt) + geom_violin(aes(x = variable, y = value)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("figures/meth_violins.png")

#density meth_diffs

meth_vals <- meth_vals[,-13]

#Cecum
cecum <- meth_vals[,metadata$segment == "cecum"]
meta_cecum <- metadata[metadata$segment == "cecum",]

diffs_cecum <- apply(cecum, 1, 
                     function(w){
                       mean(w[meta_cecum$age_group == "old"]) - 
                         mean(w[meta_cecum$age_group == "young"])})


#sigmoid
sigmoid <- meth_vals[,metadata$segment == "sigmoid"]
meta_sigmoid <- metadata[metadata$segment == "sigmoid",]

diffs_sigmoid <- apply(sigmoid, 1, 
                     function(w){
                       mean(w[meta_sigmoid$age_group == "old"]) - 
                         mean(w[meta_sigmoid$age_group == "young"])})

d2 <- data.frame(beta = c(diffs_cecum, diffs_sigmoid),
                  section = c(rep("cecum",length(diffs_cecum)), 
                              rep("sigmoid",length(diffs_sigmoid))))

ggplot(d2) + geom_density(aes(x = beta, color = section)) +
  theme_bw()
ggsave("figures/hist_betas_persection.png")


#### Run dmrseq ####

#All samples: Use segment as a covariate

set.seed(1234)
DMRsage <- dmrseq(bs=bismarkBSseq,
                  testCovariate="age_group", 
                  cutoff = 0.05, 
                  BPPARAM = MulticoreParam(3),
                  adjustCovariate = "segment",
                  maxPerms = 20,
                  maxGap = 100,
                  maxGapSmooth = 1000,
                  minNumRegion = 3)

DMRsage$beta <- -1 * DMRsage$beta #24,162
DMRsage_annot <- get_table_with_annots(DMRsage)
save(DMRsage_annot,file = "data/rdata/DMRs_age_final.RData")

df <- as.data.frame(DMRsage_annot)
write.table(df, file = "data/tables/DMRs_age_final.txt", quote = FALSE, row.names = FALSE)


# Segment: old and young together
set.seed(1234)

DMRsseg <- dmrseq(bs=bismarkBSseq,
                    testCovariate="segment",
                    adjustCovariate = "age_group",
                    #adjustCovariate = "patient",
                    cutoff = 0.05, 
                    BPPARAM = MulticoreParam(3),
                    maxPerms = 20,
                    maxGap = 100,
                    maxGapSmooth = 1000,
                    minNumRegion = 3)

DMRsage_annot <- get_table_with_annots(DMRsseg)
DMRsage_annot$beta <- -1 * DMRsage_annot$beta

save(DMRsage_annot,file = "data/rdata/DMRs_segm.RData")

df <- as.data.frame(DMRsage_annot)
write.table(df, file = "data/tables/DMRs_segm.txt", quote = FALSE, row.names = FALSE)
