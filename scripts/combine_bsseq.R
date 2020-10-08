################################################################################
# R script to combine bsseq objects from run_dmrseq_*.R scripts.
#
# dec 16 2019
################################################################################


suppressPackageStartupMessages({
  library(dmrseq)
  library(ggplot2)
})

load("data/rdata/bsseq_age_filt.RData")
bsseq_age <- bismarkBSseq
load("data/rdata/bsseq_lesions_filt.RData")

bsList <- list(bsseq_age, bismarkBSseq)
bsCombined <- combineList(bsList)

#Filter
loci.idx <- DelayedMatrixStats::rowSums2(getCoverage(bsCombined, type="Cov") >= 10 ) >= 38
bsCombined <- bsCombined[loci.idx,] #2,407,522 loci

save(bsCombined, file = "data/rdata/bsseqCombined.RData")

#save tab matrix in tab delimited format for ArrayExpress submission
# cov <- getCoverage(bsCombined, type = "Cov")
# meth <- getCoverage(bsCombined, type = "M")
# 
# write.table(cov, file = "total_reads.txt")
# write.table(meth, file = "methylated_reads.txt")
# colData(bsCombined) #have to change colnames for names in table 1 of paper
