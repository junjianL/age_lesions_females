## Explore effect sizes
# July 30 2020

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(bsseq)
  library(ggplot2)
  library(annotatr)
})


#### load combined object ####
load("data/bsseqCombined.RData")

#Get meth table
gr <- rowRanges(bsCombined)
cov <- getCoverage(bsCombined, type = "Cov")
meth <- getCoverage(bsCombined, type = "M")
seqlevels(gr) <- paste0("chr",seqlevels(gr))
meth_vals <- meth / cov
idx <- rowSums(meth_vals, na.rm=TRUE) != 0 #?
meth_vals <- meth_vals[idx,]
gr <- gr[idx]

#get regions to subset
annotsgene <- c("hg19_cpg_islands")
annotations_genes = build_annotations(genome = 'hg19', annotations = annotsgene)

grab_sites <- function(regions, gr, df){
  hits <- subjectHits(findOverlaps(regions, gr))
  dfsub <- df[hits,]
}

methvalssub <- grab_sites(annotations_genes, gr, meth_vals)

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
  scale_color_brewer(palette = "Set1") +
  geom_line(stat="density") + 
  theme_bw() +
  ggtitle("Methylation values in CpG Islands")

