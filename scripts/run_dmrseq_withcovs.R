
# dmrseq run with cov files
# oct 24 2019

suppressPackageStartupMessages({
  library(dmrseq)
  library(dplyr)
  library(ggplot2)
})


metafile <- "data/metadata.txt"
#print(metafile)

#Read in metadata
metadata <- read.delim(metafile, header = TRUE, as.is = TRUE, sep = "\t")
infile <- sprintf("data/methextract/%s_pe.CpG_report.txt", metadata$names)
infile2 <- sprintf("data/methextract/%s_pe.bismark.cov.gz", metadata$names)

#metadata$patient <- gsub("[0-9]+\\.A-", "", metadata$names)
#metadata$patient <- gsub("_[a-z]+","",metadata$patient)

metadata$age_group <- ifelse(metadata$age > 40, "old", "young")

#Read in cov files
bismarkBSseq <- read.bismark(
  files= infile2,
  rmZeroCov = TRUE,
  strandCollapse = TRUE,
  verbose = TRUE, 
  colData = metadata,
  BPPARAM = MulticoreParam(1),
  nThread = 1,
  BACKEND = "HDF5Array",
  dir = "data/HDF5_cov",
  replace = TRUE) #47,694,038

#Filter coverage low
loci.idx <- DelayedMatrixStats::rowSums2(getCoverage(bismarkBSseq, type="Cov") >= 10 ) >= 10
bismarkBSseq <- bismarkBSseq[loci.idx,] #2,798,959

#Filter coverage high
loci.idx <- DelayedMatrixStats::rowSums2(getCoverage(bismarkBSseq, type="Cov") <= 200 ) == 12
bismarkBSseq <- bismarkBSseq[loci.idx,] #2,795,901

#Filter chrom X,Y and weird ones
loci.idx <- seqnames(bismarkBSseq) %in% c(1:22)
bismarkBSseq <- bismarkBSseq[loci.idx,] #2,728,049

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

sapply(colnames(meth_vals)[1], make_bigwigs, 
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
ggsave("figures/covVsmeth.png")

#MD  

#MV
# ggplot(d) + geom_point(aes(mean_meth, var_meth), alpha = 0.2) +
#   theme_bw() +
# ggsave("figures/meanVsvar.png")

#meth per sample
pos <- paste0(seqnames(bismarkBSseq),".", start(bismarkBSseq))
meth_vals <- as.data.frame(meth_vals)
colnames(meth_vals) <- metadata$names
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

ggplot(d2) + geom_density(aes(x = beta, fill= section, 
                              color = section), alpha = 0.2) +
  theme_bw()
ggsave("figures/hist_betas_persection.png")

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

#Run dmrseq
bsc <- bismarkBSseq[,metadata$segment == "cecum"]
DMRsage_cecum <- dmrseq(bs=bsc,
                  testCovariate="age_group", 
                  cutoff = 0.05, 
                  BPPARAM = MulticoreParam(1),
                  #adjustCovariate = "patient",
                  maxGap = 100,
                  maxGapSmooth = 1000,
                  minNumRegion = 3)

save(DMRsage_cecum,file = "data/DMRs_age_cecum.RData")
load("data/DMRs_age_cecum.RData")

bss <- bismarkBSseq[,metadata$segment == "sigmoid"]
DMRsage_sig <- dmrseq(bs=bss,
                        testCovariate="age_group", 
                        cutoff = 0.05, 
                        BPPARAM = MulticoreParam(1),
                        #adjustCovariate = "patient",
                        maxGap = 100,
                        maxGapSmooth = 1000,
                        minNumRegion = 3)

save(DMRsage_cecum,file = "data/DMRs_age_sig.RData")


## DMR plots
df <- as.data.frame(mcols(DMRsage_cecum))

ggplot(df) + geom_histogram(aes(beta), bins = 50) + theme_bw()
ggsave("figures/DMRscecum_betas.png")

#Use only segment as a covariate (?)
DMRsage <- dmrseq(bs=bismarkBSseq,
                      testCovariate="age_group", 
                      cutoff = 0.05, 
                      BPPARAM = MulticoreParam(1),
                      adjustCovariate = "segment",
                      maxGap = 100,
                      maxGapSmooth = 1000,
                      minNumRegion = 3)



