#########################################################
## Plot Figure 2A and supp Fig 4A
# July 30 2020
#########################################################

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(bsseq)
  library(ggplot2)
  library(annotatr)
  library(RColorBrewer)
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

#proms
# annotsgene <- c("hg19_genes_promoters")
# annotations_genes = build_annotations(genome = 'hg19', annotations = annotsgene)
# annotations_genes <- unique(annotations_genes) #unique transcripts, not genes
# annotations_genes <- annotations_genes[!duplicated(annotations_genes$gene_id)]

grab_sites <- function(regions, gr, df){
  hits <- subjectHits(findOverlaps(regions, gr))
  dfsub <- df[hits,]
}

methvalssub <- grab_sites(annotations_genes, gr, meth_vals) #1,060,602 islands
# 400,905 proms

#### draw densities ####
bsCombined$lesion_2 <- gsub("sigmoid_|cecum_","", bsCombined$lesion)
bsCombined$lesion_2 <- factor(bsCombined$lesion_2, 
                              levels =c("Adenoma", "SSA", "Normal_Adenoma", 
                              "Normal_SSA", "old", "young"))

df <- data.frame(t(methvalssub), Tissue = bsCombined$lesion_2, Sample = bsCombined$names)

df <- reshape2::melt(df, id.vars=c("Sample", "Tissue"),
                     value.name="Methylation")

# Red: Adenoma  
# Blue: cecum_old  
# Green: cecum_young  
# Purple: Normal_Adenoma  
# Orange: Normal_SSA  
# Yellow: sigmoid_old     
# Brown: sigmoid_young   
# Pink: SSA

myColor <- c(brewer.pal(8, "Set1")[c(1,8,4,5)],"black", brewer.pal(8, "Set1")[3])
ggplot(df) + 
  scale_y_log10() +
  scale_color_manual(values=myColor) +
  geom_line(aes_string(x="Methylation", col="Tissue"),
            stat="density", adjust = 6, size = 1) + 
  theme_bw() +
  ggtitle("CpG site methylation values in CpG islands")
#ggsave("density_islands_yscalelog10.pdf")

