## Explore effect sizes
# July 16 2020

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(bsseq)
  library(ggplot2)
  library(cowplot)
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

#### Plot effect sizes ####
# in dmrseq this is defines as:
# The estimate of the mean methylation proportion is for loci i in condition s 
# is taken to be the sum of methylated reads from all samples in that condition 
# divided by the sum of all reads (i.e. the coverage) from all samples in condition s.
# This leads to the following estimate of methylation proportion difference:
# beta = prop1 - prop2 (no weighting)

get_change <- function(idx1, idx2){
prop1 <- rowSums(meth[,idx1]) / rowSums(cov[,idx1])
prop2 <- rowSums(meth[,idx2]) / rowSums(cov[,idx2])
diff <- prop1 - prop2
}

ssa <- get_change(colData(bsCombined)$lesion == "SSA", 
                  colData(bsCombined)$lesion == "Normal_SSA")

cadn <- get_change(colData(bsCombined)$lesion == "Adenoma",
  colData(bsCombined)$lesion == "Normal_Adenoma")

age <- get_change(colData(bsCombined)$lesion %in% c("cecum_old", "sigmoid_old"), 
  colData(bsCombined)$lesion %in% c("cecum_young", "sigmoid_young"))

seg <- get_change(colData(bsCombined)$lesion %in% c("cecum_old", "cecum_young"), 
  colData(bsCombined)$lesion %in% c("sigmoid_old", "sigmoid_young"))

df <- data.frame("SSA/P" = ssa,
                 cADN = cadn,
                 Age = age,
                 Segment = seg)

# pairs

custom_func <- function(data, mapping){
  ggplot(data = data, mapping = mapping) +
    geom_bin2d() +
    scale_fill_distiller(palette='RdBu', trans='log10') +
    scale_x_continuous(limits = c(-0.6,0.8)) +
    scale_y_continuous(limits = c(-0.6,0.9)) 
}

custom_func2 <- function(data, mapping){
  ggplot(data = data, mapping = mapping) +
    geom_density() +
    scale_x_continuous(limits = c(-0.6,0.8)) +
    scale_y_sqrt(limits = c(0,120))
}

GGally::ggpairs(df, 
                upper = list(continuous = custom_func2),
                lower = list(continuous = custom_func),
                diag = list(na = "naDiag"))

## within genomic compartments

#promoters
annotsgene <- c("hg19_genes_promoters")
annotations_genes = build_annotations(genome = 'hg19', annotations = annotsgene)
annotations_genes <- unique(annotations_genes) #unique transcripts, not genes
annotations_genes <- annotations_genes[!duplicated(annotations_genes$gene_id)]

#islands
annotsgene <- c("hg19_cpg_islands")
annotations_genes <- build_annotations(genome = 'hg19', annotations = annotsgene)

grab_sites <- function(regions, gr, df){
hits <- subjectHits(findOverlaps(regions, gr))
dfsub <- df[hits,]
#
}

df_proms <- grab_sites(annotations_genes, gr, df) #1,064,915
df_proms <- df_proms[rowSums(df_proms) > 0,] #653,484

GGally::ggpairs(df_proms,
                upper = NULL,
                lower = list(continuous = custom_func),
                diag = list(continuous = custom_func2)) +
  theme_bw()


#### Plot meth values ####

#choose DMR lists
cutoff <- 0.05
load("data/rdata/DMRs_lesions_SSA.RData")
ssa <- DMRsles_annot[DMRsles_annot$qval <= cutoff & DMRsles_annot$beta > 0]

load("data/rdata/DMRs_lesions_cADN.RData")
cadn <- DMRsles_annot[DMRsles_annot$qval <= cutoff & DMRsles_annot$beta > 0]

load("data/rdata/DMRs_age_final.RData")
age <- DMRsage_annot[DMRsage_annot$qval <= cutoff & DMRsage_annot$beta > 0]

load("data/rdata/DMRs_segm.RData")
seg <- DMRsage_annot[DMRsage_annot$qval <= cutoff & DMRsage_annot$beta < 0]

# get df for each comparison and draw plot
plot_meth <- function(regs, idx1, idx2, labx, laby){
  methsub <- grab_sites(regs, gr, meth_vals)
  #covsub <- grab_sites(regs, gr, cov)
  
  #prop1 <- rowSums(methsub[,idx1]) / rowSums(covsub[,idx1])
  prop1 <- rowMeans(methsub[,idx1], na.rm = TRUE)
  prop2 <- rowMeans(methsub[,idx2], na.rm = TRUE)
  dat <- data.frame(prop1, prop2)
  
  print(mean(prop1))
  print(mean(prop2))
  
  #return(dat)
  p <- ggplot(data = dat, aes(prop1, prop2)) +
    geom_bin2d() +
    #geom_density_2d(color = "black") +
    scale_fill_distiller(palette='RdBu', trans='log10', limits = c(1,10605)) +
    geom_abline() +
    ylab(laby) + xlab(labx) +
    theme_bw()
  p
  
}


a <- plot_meth(ssa, bsCombined$lesion == "Normal_SSA",
               bsCombined$lesion == "SSA",
               "Normal SSA", "SSA/P")


b <- plot_meth(cadn, colData(bsCombined)$lesion == "Normal_Adenoma",
               colData(bsCombined)$lesion == "Adenoma",
               "Normal cADN", "cADN")

c <- plot_meth(age, 
               colData(bsCombined)$lesion %in% c("cecum_young", "sigmoid_young"),
               colData(bsCombined)$lesion %in% c("cecum_old", "sigmoid_old"), 
               "Young", "Old")
d <- plot_meth(seg, 
               colData(bsCombined)$lesion %in% c("cecum_young", "cecum_old"),
               colData(bsCombined)$lesion %in% c("sigmoid_young", "sigmoid_old"),
               "Cecum","Sigmoid")


cowplot::plot_grid(a,b,c,d, labels = "AUTO", ncol = 2)
ggsave("scatter_methvals_inrespectice_hyperDMRs_2dens.pdf", width = 8, height = 6)
