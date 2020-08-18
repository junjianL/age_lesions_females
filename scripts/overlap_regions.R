## Overlap Regions
# jul 2020

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(bsseq)
})


#### Get regions that are uniquely in pre-lesions ####

cutoff <- 0.05

# get all combinations of lesion comparisons

load("data/rdata/DMRs_lesions_3.RData")
full <- sort(DMRsles_annot[DMRsles_annot$qval <= cutoff & DMRsles_annot$beta > 0])
full$comparison <- "SSA+cADN"
load("data/rdata/DMRs_lesions_SSA.RData")
ssa <- sort(DMRsles_annot[DMRsles_annot$qval <= cutoff & DMRsles_annot$beta > 0])
ssa$comparison <- "SSA"
load("data/rdata/DMRs_lesions_cADN.RData")
cadn <- sort(DMRsles_annot[DMRsles_annot$qval <= cutoff & DMRsles_annot$beta > 0])
cadn$comparison <- "cADN"
rm(DMRsles_annot)

les_source <- c(full,ssa,cadn)
o <- order(les_source$qval)
les_source <- les_source[o]
les <- reduce(les_source, with.revmap=TRUE)

# get all combinations of age comparisons

load("data/rdata/DMRs_age_final.RData")
full <- sort(DMRsage_annot[DMRsage_annot$qval <= cutoff])
load("data/rdata/DMRs_age_sig_2.RData")
sig <- sort(DMRsage_sig_annot[DMRsage_sig_annot$qval <= cutoff])
load("data/rdata/DMRs_age_cecum_2.RData")
cec <- sort(DMRsage_cecum_annot[DMRsage_cecum_annot$qval <= cutoff])

age <- c(full, sig, cec)
age <- reduce(age)

#get all combinations of segment comparisons

load("data/rdata/DMRs_segm.RData")
full <- sort(DMRsage_annot[DMRsage_annot$qval <= cutoff])
load("data/rdata/DMRs_young.RData")
young <- sort(DMRsage_annot[DMRsage_annot$qval <= cutoff])
load("data/rdata/DMRs_old.RData")
old <- sort(DMRsage_annot[DMRsage_annot$qval <= cutoff])
rm(DMRsage_annot, DMRsage_cecum_annot, DMRsage_sig_annot)

seg <- c(full, young, old)
seg <- reduce(seg)

age_seg <- c(age, seg)
age_seg <- reduce(age_seg)

#get lesion regions without age or segment signal

unique_less <- subsetByOverlaps(les_source, age_seg, invert = TRUE)
unique_less
save(unique_less,file = sprintf("data/rdata/uniqueDMRs_lesion_%g.RData", cutoff))

df <- as.data.frame(unique_less)
write.table(df, file = sprintf("data/tables/uniqueDMRs_lesion_%g.txt",cutoff), 
            quote = FALSE, row.names = FALSE)


# merged
unique_less <- subsetByOverlaps(les, age_seg, invert = TRUE)
unique_less
save(unique_less,file = sprintf("data/rdata/uniqueDMRs_lesion_merged_%g.RData", cutoff))

df <- as.data.frame(unique_less)
write.table(df, file = sprintf("data/tables/uniqueDMRs_lesion_merged_%g.txt",cutoff), 
            quote = FALSE, row.names = FALSE)


#### filter regions ####

load("data/rdata/bsseq_lesions_filt.RData")
#load("data/rdata/uniqueDMRs_lesion_merged_0.05.RData")

#Get meth table
gr <- rowRanges(bismarkBSseq)
cov <- getCoverage(bismarkBSseq, type = "Cov")
meth <- getCoverage(bismarkBSseq, type = "M")
seqlevels(gr) <- paste0("chr",seqlevels(gr))


#summarize meth val per region
library(plyranges)
library(dplyr)
mcols(gr) <- meth
hits <- findOverlaps(unique_less, gr)
gr$DMR <- NA
gr[subjectHits(hits)]$DMR <- queryHits(hits)
gr_sub <- gr[!is.na(gr$DMR)]

#methylation matrix
meth_dmr <- gr_sub %>% 
  group_by(DMR) %>% 
  summarise_at(
    colnames(meth), mean, na.rm=TRUE
  ) %>% 
  as.matrix()

#cov matrix
mcols(gr) <- cov
gr$DMR <- NA
gr[subjectHits(hits)]$DMR <- queryHits(hits)
gr_sub <- gr[!is.na(gr$DMR)]

cov_dmr <- gr_sub %>% 
  group_by(DMR) %>% 
  summarise_at(
    colnames(cov), mean, na.rm=TRUE
  ) %>% 
  as.matrix()

# 1. Filter by baseline level

meth_val <- meth_dmr / cov_dmr
idxnorm <- grepl("Normal_",colData(bismarkBSseq)$lesion)
subnorm <- as.matrix(meth_val[,idxnorm])
idx <- rowMeans(subnorm, na.rm = TRUE) < 0.1    

# 2. Filter by effect size
get_change <- function(idx1, idx2){
  prop1 <- rowSums(meth_dmr[,idx1]) / rowSums(cov_dmr[,idx1])
  prop2 <- rowSums(meth_dmr[,idx2]) / rowSums(cov_dmr[,idx2])
  diff <- prop1 - prop2
}

lesions_diff <- get_change(!idxnorm, idxnorm)
ssa <- get_change(colData(bismarkBSseq)$lesion == "SSA", 
                  colData(bismarkBSseq)$lesion == "Normal_SSA")

cadn <- get_change(colData(bismarkBSseq)$lesion == "Adenoma",
                   colData(bismarkBSseq)$lesion == "Normal_Adenoma")
  
cutoff <- 0.8
idx2 <- lesions_diff > cutoff | ssa > cutoff | cadn > cutoff

# Do filter
sub_unique <- unique_less[idx & idx2]
save(sub_unique, file = "data/rdata/unique_lesions_filt.RData")

dfgr <- as.data.frame(sub_unique)
write.csv(dfgr, file = "data/tables/uniqueDMRs_lesion_merged_filt.txt", 
          quote = FALSE, row.names = FALSE)
