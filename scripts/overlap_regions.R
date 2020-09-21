## Overlap Regions
# jul 2020

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(bsseq)
  library(plyranges)
  library(dplyr)
  library(annotatr)
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

les$qvals <- lapply(les$revmap, function(u){
  les_source$qval[u] 
})

les$comparison <- lapply(les$revmap, function(u){
  les_source$comparison[u] 
})

# get age comparison

load("data/rdata/DMRs_age_final.RData")
age <- sort(DMRsage_annot[DMRsage_annot$qval <= cutoff])

#get segment comparison

load("data/rdata/DMRs_segm.RData")
seg <- sort(DMRsage_annot[DMRsage_annot$qval <= cutoff])

age_seg <- c(age, seg)
age_seg <- reduce(age_seg)

#get lesion regions without age or segment signal
# merged
unique_less <- subsetByOverlaps(les, age_seg, invert = TRUE)
unique_less

#get unique age regions
les_seg <- c(les, seg)
les_seg <- reduce(les_seg)

unique_age <- subsetByOverlaps(age, les_seg, invert = TRUE)
df <- as.data.frame(unique_age)
write.table(df, file = "data/tables/uniqueDMRs_age.txt",
            quote = FALSE, row.names = FALSE)

#get unique segment regions
les_age <- c(les, age)
les_age <- reduce(les_age)

unique_seg <- subsetByOverlaps(seg, les_age, invert = TRUE)
df <- as.data.frame(unique_seg)
write.table(df, file = "data/tables/uniqueDMRs_segment.txt",
            quote = FALSE, row.names = FALSE)


#### filter regions and get tumor-unique regions ####

load("data/rdata/bsseq_lesions_filt.RData")

#Get meth table
gr <- rowRanges(bismarkBSseq)
cov <- getCoverage(bismarkBSseq, type = "Cov")
meth <- getCoverage(bismarkBSseq, type = "M")
seqlevels(gr) <- paste0("chr",seqlevels(gr))

#summarize meth val per region
mcols(gr) <- meth
hits <- findOverlaps(unique_less, gr)
gr$DMR <- NA
gr[subjectHits(hits)]$DMR <- queryHits(hits)
gr_sub <- gr[!is.na(gr$DMR)]

#meth cov matrix
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
#plot(lesions_diff, lesions_diff2)

ssa <- get_change(colData(bismarkBSseq)$lesion == "SSA", 
                  colData(bismarkBSseq)$lesion == "Normal_SSA")

cadn <- get_change(colData(bismarkBSseq)$lesion == "Adenoma",
                   colData(bismarkBSseq)$lesion == "Normal_Adenoma")
  
cutoff <- 0.8
idx2 <- lesions_diff > cutoff | ssa > cutoff | cadn > cutoff

# Do filter
sub_unique <- unique_less[idx & idx2]

#add original qvalues to matrix
qvalmat <- data.frame(matrix(NA, nrow = length(sub_unique), ncol=3))
colnames(qvalmat) <- c("SSA+cADN", "SSA", "cADN")
head(qvalmat)  

for (u in 1:length(sub_unique)){
  comps <- sub_unique$comparison[[u]]
  comps <- comps[!duplicated(comps)]
  
  ord <- match(comps, colnames(qvalmat))
  qvalmat[u,ord] <- sub_unique$qvals[[u]][!duplicated(comps)]
}

head(qvalmat)  
mcols(sub_unique) <- qvalmat

#annotate
source("scripts/get_table_with_annots.R")
sub_uniqueannot <- get_table_with_annots(sub_unique, suffix = "chr")

#save object and file
save(sub_uniqueannot, file = "data/rdata/unique_lesions_filt.RData")

dfgr <- as.data.frame(sub_uniqueannot)
head(dfgr)
write.table(dfgr, file = "data/tables/uniqueDMRs_lesion_merged_filt.txt", 
          quote = FALSE, row.names = FALSE, sep = "\t")
