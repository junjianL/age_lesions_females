
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
metafile <- "data/metadata.txt"
#print(metafile)

#Read in metadata
metadata <- read.delim(metafile, header = TRUE, as.is = TRUE, sep = "\t")
infile <- sprintf("data/methextract/%s_pe.CpG_report.txt", metadata$names)
#infile2 <- sprintf("data/methextract/%s_pe.bismark.cov.gz", metadata$names)

#metadata$patient <- gsub("[0-9]+\\.A-", "", metadata$names)
#metadata$patient <- gsub("_[a-z]+","",metadata$patient)

metadata$age_group <- ifelse(metadata$age > 50, "old", "young")


#Read in cov files
bismarkBSseq <- read.bismark(files = infile,
                             #files= infile2,
                             rmZeroCov = TRUE,
                             strandCollapse = TRUE,
                             verbose = TRUE, 
                             colData = metadata,
                             BPPARAM = MulticoreParam(1),
                             nThread = 1,
                             BACKEND = "HDF5Array",
                             dir = "data/HDF5_cov",
                             replace = TRUE) #26,660,520

save(bismarkBSseq, file = "data/bsseq_cpgreport.RData")
#load("data/bsseq_cpgreport.RData")

#Filter coverage low
loci.idx <- DelayedMatrixStats::rowSums2(getCoverage(bismarkBSseq, type="Cov") >= 10 ) >= 10
bismarkBSseq <- bismarkBSseq[loci.idx,] #2,481,876

#Filter chrom X,Y and weird ones
loci.idx <- seqnames(bismarkBSseq) %in% c(1:22)
bismarkBSseq <- bismarkBSseq[loci.idx,] #2,418,647

#Filter/reduce coverage high (test)
#loci.idx <- DelayedMatrixStats::rowSums2(getCoverage(bismarkBSseq, type="Cov") <= 200 ) == 36
#bismarkBSseq <- bismarkBSseq[loci.idx,] 

#or from BiSeq
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
       chromsizesFile = "/home/Shared_taupo/steph/reference/hg19.chrom.sizes.modGL")


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

#histogram meth_diffs

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


#cov per sample
# cov_vals_melt <- reshape2::melt(cov, measure.vars = metadata$names)
# ggplot(cov_vals_melt) + geom_violin(aes()) +
#   theme_bw()
# ggsave("figures/cov_violins.png")


#### Run dmrseq ####
#scale beta func

scale_beta <- function(u){
  min_u <- min(abs(u))
  max_u <- max(abs(u))
  sign(u) * (abs(u)-min_u)/(max_u-min_u)
  # min_u <- min(u)
  # max_u <- max(u)
  # (u-min_u)/(max_u-min_u)
}

#Cecum

colData(bismarkBSseq)$age_group <- relevel(factor(metadata$age_group), ref = "old")
bsc <- bismarkBSseq[,metadata$segment == "cecum"]
DMRsage_cecum <- dmrseq(bs=bsc,
                        testCovariate="age_group", 
                        cutoff = 0.05, 
                        BPPARAM = MulticoreParam(1),
                        #adjustCovariate = "patient",
                        maxGap = 100,
                        maxGapSmooth = 1000,
                        minNumRegion = 3)
#DMRsage_cecum <- DMRsage_cecum[DMRsage_cecum$qval < 0.05] #697
DMRsage_cecum$beta <- -1 * DMRsage_cecum$beta #31,145
#DMRsage_cecum$meth <- scale_beta(DMRsage_cecum$beta)
DMRsage_cecum_annot <- get_table_with_annots(DMRsage_cecum)

hist(DMRsage_cecum$beta)
hist(DMRsage_cecum$meth)

save(DMRsage_cecum_annot,file = "data/DMRs_age_cecum_2.RData")

df <- as.data.frame(DMRsage_cecum_annot)
write.table(df, file = "data/DMRs_age_cecum_2.txt", quote = FALSE, row.names = FALSE)

#Sigmoid

bss <- bismarkBSseq[,metadata$segment == "sigmoid"]
DMRsage_sig <- dmrseq(bs=bss,
                      testCovariate="age_group", 
                      cutoff = 0.05, 
                      BPPARAM = MulticoreParam(1),
                      #adjustCovariate = "patient",
                      maxGap = 100,
                      maxGapSmooth = 1000,
                      minNumRegion = 3)
#DMRsage_sig <- DMRsage_sig[DMRsage_sig$qval < 0.05] #112
DMRsage_sig$beta <- -1 * DMRsage_sig$beta #37,382
DMRsage_sig$meth <- scale_beta(DMRsage_sig$beta)
DMRsage_sig_annot <- get_table_with_annots(DMRsage_sig)

save(DMRsage_sig_annot,file = "data/DMRs_age_sig_2.RData")

df <- as.data.frame(DMRsage_sig_annot)
write.table(df, file = "data/DMRs_age_sigmoid_2.txt", quote = FALSE, row.names = FALSE)

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
DMRsage$meth <- scale_beta(DMRsage$beta)
DMRsage_annot <- get_table_with_annots(DMRsage)
save(DMRsage_annot,file = "data/DMRs_age_final.RData")

df <- as.data.frame(DMRsage_annot)
write.table(df, file = "data/DMRs_age_final.txt", quote = FALSE, row.names = FALSE)

# Only old 
set.seed(1234)
bso <- bismarkBSseq[,metadata$age_group == "old"]
DMRsage_o <- dmrseq(bs=bso,
                  testCovariate="segment", 
                  cutoff = 0.05, 
                  BPPARAM = MulticoreParam(3),
                  maxPerms = 20,
                  maxGap = 100,
                  maxGapSmooth = 1000,
                  minNumRegion = 3)

DMRsage_annot <- get_table_with_annots(DMRsage_o)
save(DMRsage_annot,file = "data/DMRs_old.RData")

df <- as.data.frame(DMRsage_annot)
write.table(df, file = "data/DMRs_old.txt", quote = FALSE, row.names = FALSE)

# Only young
set.seed(1234)
bsy <- bismarkBSseq[,metadata$age_group == "young"]
DMRsage_y <- dmrseq(bs=bsy,
                    testCovariate="segment", 
                    cutoff = 0.05, 
                    BPPARAM = MulticoreParam(3),
                    maxPerms = 20,
                    maxGap = 100,
                    maxGapSmooth = 1000,
                    minNumRegion = 3)

DMRsage_annot <- get_table_with_annots(DMRsage_y)
save(DMRsage_annot,file = "data/DMRs_young.RData")

df <- as.data.frame(DMRsage_annot)
write.table(df, file = "data/DMRs_young.txt", quote = FALSE, row.names = FALSE)

# old and young together
set.seed(1234)

DMRsseg <- dmrseq(bs=bismarkBSseq,
                    testCovariate="segment",
                    adjustCovariate = "age_group",
                    cutoff = 0.05, 
                    BPPARAM = MulticoreParam(3),
                    maxPerms = 20,
                    maxGap = 100,
                    maxGapSmooth = 1000,
                    minNumRegion = 3)

DMRsage_annot <- get_table_with_annots(DMRsseg)
save(DMRsage_annot,file = "data/DMRs_segm.RData")

df <- as.data.frame(DMRsage_annot)
write.table(df, file = "data/DMRs_segm.txt", quote = FALSE, row.names = FALSE)
