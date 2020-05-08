
################################################################################
# R script to read-in CpG-report files from bismark, limit coverage, and 
# generate bigwig files
#
# Stephany Orjuela, April 2020
################################################################################

args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

print(files)
print(chromSizes)

suppressPackageStartupMessages({
  library(dmrseq)
})

#infile <- sprintf("data/methextract_lesions/%s_pe.CpG_report.txt", metadata$names)


#Read in cov files
bismarkBSseq <- read.bismark(files = files,
                             rmZeroCov = TRUE,
                             strandCollapse = TRUE,
                             verbose = TRUE, 
                             #colData = metadata,
                             BPPARAM = MulticoreParam(1),
                             nThread = 1,
                             #BACKEND = "HDF5Array",
                             #dir = "data/HDF5_lesion",
                             replace = TRUE) 


#Filter coverage low
loci.idx <- DelayedMatrixStats::rowSums2(getCoverage(bismarkBSseq, type="Cov") >= 10 ) >= 1
bismarkBSseq <- bismarkBSseq[loci.idx,] 

#Filter chrom X,Y and weird ones
loci.idx <- seqnames(bismarkBSseq) %in% c(1:22)
bismarkBSseq <- bismarkBSseq[loci.idx,] 

#Filter/reduce coverage high (code from BiSeq)

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
                   stringsAsFactors = FALSE,
                   row.names = NULL)
  
  #remove NAs 
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
       chromsizesFile = chromSizes)