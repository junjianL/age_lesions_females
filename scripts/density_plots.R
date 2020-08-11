## Explore effect sizes
# July 30 2020

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(bsseq)
  library(ggplot2)
  library(annotatr)
})


#### load combined object ####
load("data/rdata/bsseqCombined.RData")

#Get meth table
gr <- rowRanges(bsCombined)
cov <- getCoverage(bsCombined, type = "Cov")
meth <- getCoverage(bsCombined, type = "M")
seqlevels(gr) <- paste0("chr",seqlevels(gr))
meth_vals <- meth / cov

#remove empty rows in all samples
idx <- rowSums(meth_vals, na.rm=TRUE) != 0 #?
meth_vals <- meth_vals[idx,]
gr <- gr[idx] #2,402,537

#get regions to subset
annotsgene <- c("hg19_cpg_islands")
annotations_genes = build_annotations(genome = 'hg19', annotations = annotsgene)

grab_sites <- function(regions, gr, df){
  hits <- subjectHits(findOverlaps(regions, gr))
  dfsub <- df[hits,]
}

methvalssub <- grab_sites(annotations_genes, gr, meth_vals) #1,060,602      48

#### draw densities ####
bsCombined$lesion_2 <- gsub("sigmoid_|cecum_","", bsCombined$lesion)
bsCombined$lesion_2 <- factor(bsCombined$lesion_2, 
                              levels =c("Adenoma", "SSA", "Normal_Adenoma", 
                              "Normal_SSA", "old", "young"))

df <- data.frame(t(methvalssub), Tissue = bsCombined$lesion_2, Sample = bsCombined$names)

df <- reshape2::melt(df, id.vars=c("Sample", "Tissue"),
                     value.name="Methylation")

ggplot(df, aes_string(x="Methylation", 
                      col="Tissue"), 
       fill=NULL) + 
  scale_x_sqrt() +
  scale_y_sqrt() +
  #expand_limits(x = 0, y = 0 ) + 
  scale_color_brewer(palette = "Set1") +
  geom_line(stat="density") + 
  theme_bw() +
  ggtitle("CpG site methylation values in CpG Islands")


### use lesion DMRs
load("data/rdata/DMRs_lesions_3.RData")
load("data/rdata/DMRs_age_final.RData")

methvalsdmr <- grab_sites(DMRsage_annot, gr, meth_vals) 
methvalsdmr <- grab_sites(DMRsles_annot, gr, meth_vals)

df <- data.frame(t(methvalsdmr), Tissue = bsCombined$lesion_2, Sample = bsCombined$names)

df <- reshape2::melt(df, id.vars=c("Sample", "Tissue"),
                     value.name="Methylation")

ggplot(df, aes_string(x="Methylation", 
                      col="Tissue"), 
       fill=NULL) + 
  scale_x_sqrt() +
  scale_y_sqrt() +
  #expand_limits(x = 0, y = 0 ) + 
  scale_color_brewer(palette = "Set1") +
  geom_line(stat="density") + 
  theme_bw() +
  ggtitle("CpG site methylation values in lesion DMRs")
  #ggtitle("CpG site methylation values in age DMRs")