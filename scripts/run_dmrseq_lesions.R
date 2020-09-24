## Run dmrseq lesions
# nov 25 2019

suppressPackageStartupMessages({
  library(dmrseq)
  library(dplyr)
  library(ggplot2)
  library(BiocParallel)
})

source("scripts/get_table_with_annots.R")
metafile <- "metadata_lesions.txt" #36 samples, 18 patients
#print(metafile)

#Read in metadata
metadata <- read.delim(metafile, header = TRUE, as.is = TRUE, sep = "\t")
infile <- sprintf("data/methextract_lesions/%s_pe.CpG_report.txt", metadata$names)

metadata$age_group <- ifelse(metadata$age > 50, "old", "young")
metadata$state <- ifelse(grepl("Normal", metadata$lesion), "Normal", "Cancer")


#Read in cov files
bismarkBSseq <- read.bismark(files = infile,
                             rmZeroCov = TRUE,
                             strandCollapse = TRUE,
                             verbose = TRUE, 
                             colData = metadata,
                             BPPARAM = MulticoreParam(1),
                             nThread = 1,
                             BACKEND = "HDF5Array",
                             dir = "data/HDF5_lesion",
                             replace = TRUE) #27,546,598

save(bismarkBSseq, file = "data/rdata/bsseq_lesions_cpgreport.RData")

#Filter coverage low
loci.idx <- DelayedMatrixStats::rowSums2(getCoverage(bismarkBSseq, type="Cov") >= 10 ) >= 28
bismarkBSseq <- bismarkBSseq[loci.idx,] #2,773,805

#Filter chrom X,Y and weird ones
loci.idx <- seqnames(bismarkBSseq) %in% c(1:22)
bismarkBSseq <- bismarkBSseq[loci.idx,] #2,703,961

#Filter/reduce coverage high

cov <- getCoverage(bismarkBSseq, type = "Cov")
meth <- getCoverage(bismarkBSseq, type = "M")
ind.cov <- cov > 0
quant <- quantile(cov[ind.cov], 0.9)
quant #92 

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

save(bismarkBSseq, file = "data/rdata/bsseq_lesions_filt.RData")

#### Generate bigwigs ####
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


#### Diagnostic plots ####
mean_meth <- rowMeans(meth_vals) 
mean_cov <- rowMeans(cov)
var_meth <- rowVars(as.matrix(meth_vals))

d <- data.frame(mean_meth = mean_meth,
                mean_cov =  mean_cov,
                #diff = diffs,
                var_meth = var_meth)

ggplot(d) + aes(x=mean_cov, y=mean_meth) + 
  geom_bin2d() + 
  scale_fill_distiller(palette='RdBu', trans='log10') + 
  scale_x_continuous(breaks=scales::pretty_breaks()) +
  theme_bw()

ggsave("figures/covVsmeth_lesions2.png")

#MV
ggplot(d) + aes(x=mean_cov, y=var_meth) + 
  geom_bin2d() + 
  scale_fill_distiller(palette='RdBu', trans='log10') + 
  scale_x_continuous(breaks=scales::pretty_breaks()) +
  theme_bw()
ggsave("figures/meanVsvar_lesions2.png")

#meth per sample
pos <- paste0(seqnames(bismarkBSseq),".", start(bismarkBSseq))
meth_vals <- as.data.frame(meth_vals)
meth_vals$pos <- pos
meth_vals_melt <- reshape2::melt(meth_vals, id.vars = "pos", measure.vars = metadata$names)

ggplot(meth_vals_melt) + geom_violin(aes(x = variable, y = value)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("figures/meth_violins_lesions.png")


#### Run dmrseq ####
set.seed(1234)
DMRsles <- dmrseq(bs=bismarkBSseq,
                  testCovariate="state", 
                  cutoff = 0.05, 
                  BPPARAM = MulticoreParam(3),
                  adjustCovariate = "patient",
                  maxPerms = 20,
                  maxGap = 100,
                  maxGapSmooth = 1000,
                  minNumRegion = 3) #55,186

#Filter
DMRsles_filt <- DMRsles[!is.na(DMRsles$qval)] #55,128
DMRsles_filt$beta <- -1 * DMRsles_filt$beta
DMRsles_filt$meth <- scale_beta(DMRsles_filt$beta)
DMRsles_annot <- get_table_with_annots(DMRsles_filt)

#save files
save(DMRsles_annot,file = "data/DMRs_lesions_3.RData")

df <- as.data.frame(DMRsles_annot)
write.table(df, file = "data/DMRs_lesions_3.txt", quote = FALSE, row.names = FALSE)

#### only SSA ####
bs.ssa <- bismarkBSseq[,grepl("SSA",colData(bismarkBSseq)$lesion)]
set.seed(1234)
DMRsles.ssa <- dmrseq(bs=bs.ssa,
                  testCovariate="lesion", 
                  cutoff = 0.05, 
                  BPPARAM = MulticoreParam(4),
                  adjustCovariate = "patient",
                  maxPerms = 20,
                  maxGap = 100,
                  maxGapSmooth = 1000,
                  minNumRegion = 3) 

#Filter
DMRsles_filt <- DMRsles.ssa[!is.na(DMRsles.ssa$qval)]
DMRsles_annot <- get_table_with_annots(DMRsles_filt)

#add number of CpGs
DMRsles_annot$numCpGs <- countOverlaps(DMRsles_filt, rowRanges(bismarkBSseq))


#save files
save(DMRsles_annot,file = "data/DMRs_lesions_SSA.RData")
df <- as.data.frame(DMRsles_annot)
write.table(df, file = "data/DMRs_lesions_SSA.txt", quote = FALSE, row.names = FALSE)

#### only cADN ####
bs.cadn <- bismarkBSseq[,grepl("Adenoma",colData(bismarkBSseq)$lesion)]
colData(bs.cadn)$lesion <- factor(colData(bs.cadn)$lesion, levels = c("Normal_Adenoma", "Adenoma")) 
set.seed(1234)
DMRsles.cadn <- dmrseq(bs=bs.cadn,
                      testCovariate="lesion", 
                      cutoff = 0.05, 
                      BPPARAM = MulticoreParam(4),
                      adjustCovariate = "patient",
                      maxPerms = 20,
                      maxGap = 100,
                      maxGapSmooth = 1000,
                      minNumRegion = 3) 

#Filter
DMRsles_filt <- DMRsles.cadn[!is.na(DMRsles.cadn$qval)]
DMRsles_annot <- get_table_with_annots(DMRsles_filt)

#save files
save(DMRsles_annot,file = "data/DMRs_lesions_cADN.RData")
df <- as.data.frame(DMRsles_annot)
write.table(df, file = "data/DMRs_lesions_cADN.txt", quote = FALSE, row.names = FALSE)

