## Run dmrseq lesions
# nov 25 2019

suppressPackageStartupMessages({
  library(dmrseq)
  library(dplyr)
  library(ggplot2)
})


metafile <- "data/metadata_serrated.txt" #36 samples, 18 patients
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

save(bismarkBSseq, file = "data/bsseq_lesions_cpgreport.RData")

#Filter coverage low
loci.idx <- DelayedMatrixStats::rowSums2(getCoverage(bismarkBSseq, type="Cov") >= 10 ) >= 28
bismarkBSseq <- bismarkBSseq[loci.idx,] #2,773,805

#Filter coverage high
loci.idx <- DelayedMatrixStats::rowSums2(getCoverage(bismarkBSseq, type="Cov") <= 200 ) == 36
bismarkBSseq <- bismarkBSseq[loci.idx,] #2,628,330

#Filter chrom X,Y and weird ones
loci.idx <- seqnames(bismarkBSseq) %in% c(1:22)
bismarkBSseq <- bismarkBSseq[loci.idx,] #2,563,041


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
ggplot(d) + geom_point(aes(x= mean_cov, y = mean_meth), alpha = 0.2) + 
  #geom_density_2d() +
  theme_bw()
ggsave("figures/covVsmeth_lesions.png")

#MD  

#MV
# ggplot(d) + geom_point(aes(mean_meth, var_meth), alpha = 0.2) +
#   theme_bw() +
# ggsave("figures/meanVsvar.png")

#meth per sample
pos <- paste0(seqnames(bismarkBSseq),".", start(bismarkBSseq))
meth_vals <- as.data.frame(meth_vals)
meth_vals$pos <- pos
meth_vals_melt <- reshape2::melt(meth_vals, id.vars = "pos", measure.vars = metadata$names)

ggplot(meth_vals_melt) + geom_violin(aes(x = variable, y = value)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("figures/meth_violins_lesions.png")

#Run dmrseq
DMRsles <- dmrseq(bs=bismarkBSseq,
                  testCovariate="state", 
                  cutoff = 0.05, 
                  BPPARAM = MulticoreParam(1),
                  adjustCovariate = "patient",
                  maxGap = 100,
                  maxGapSmooth = 1000,
                  minNumRegion = 3)
save(DMRsage_cecum,file = "data/DMRs_lesions.RData")

