## Explore effect sizes
# July 16 2020

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(bsseq)
  library(ggplot2)
})


#### load combined object ####
load("data/bsseqCombined.RData")


#Get meth table
gr <- rowRanges(bsCombined)
cov <- getCoverage(bsCombined, type = "Cov")
meth <- getCoverage(bsCombined, type = "M")
seqlevels(gr) <- paste0("chr",seqlevels(gr))

#Fix some annotations
colData(bsCombined)$state <- ifelse(is.na(colData(bsCombined)$state), 
                                    "Normal", colData(bsCombined)$state)

colData(bsCombined)$state <- ifelse(grepl("SSA|Adenoma", colData(bsCombined)$lesion), 
                                    colData(bsCombined)$lesion, colData(bsCombined)$state)

colData(bsCombined)$state <- ifelse(grepl("Adenoma", colData(bsCombined)$state), 
                                    gsub("Adenoma","cADN", colData(bsCombined)$state), 
                                    colData(bsCombined)$state)

colData(bsCombined)$state <- ifelse(grepl("SSA", colData(bsCombined)$state), 
                                    gsub("SSA","SSA/P", colData(bsCombined)$state), 
                                    colData(bsCombined)$state)

diagnosis <- ifelse(grepl("SSA|cADN", colData(bsCombined)$state), 
                    colData(bsCombined)$state, "healthy")

colData(bsCombined)$diagnosis <- factor(gsub("Normal_","",diagnosis), 
                                        levels = c("healthy", "cADN", "SSA/P"))

colData(bsCombined)$tissue <- factor(ifelse(grepl("Normal", colData(bsCombined)$state), 
                                            "normal mucosa","lesion"), levels = c("normal mucosa","lesion"))


seg <- ifelse(colData(bsCombined)$segment == "C","cecum", colData(bsCombined)$segment)
colData(bsCombined)$seg <- factor(ifelse(seg == "A","ascend", seg), levels = c("sigmoid", "ascend", "cecum"))


age <- ifelse(colData(bsCombined)$age < 40, "<=40", "41-70")
age <- ifelse(colData(bsCombined)$age > 70, ">70", age)
colData(bsCombined)$age_group <- factor(age, levels = c("<=40", "41-70", ">70"))

colnames(meth_vals) <- paste0(colData(bsCombined)$patient,".",colData(bsCombined)$lesion)
mcols(gr) <- meth_vals
seqlevels(gr) <- paste0("chr",seqlevels(gr))

#### Plot violins of effect sizes ####
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

e3 <- tidyr::gather(df, Contrast, Change)
e3$Contrast <- factor(e3$Contrast, levels = c("SSA.P", "cADN", "Age", "Segment"))



ggplot(e3) +
  geom_violin(aes(x=Contrast, y=asin(2*Change-1)) , fill = "yellow"
              ) +# scale = "width") +
  geom_boxplot(aes(x=Contrast, y=asin(2*Change-1)), color = "grey",
               alpha = 0) +
  #scale_y_continuous(trans='log10') +
  #scale_fill_manual(values = myColor) +
  #scale_color_manual(values = myColor) +
  theme_bw() +
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        text = element_text(size = 15),
        legend.position = "none")

# pairs

custom_func <- function(data, mapping){
  ggplot(data = data, mapping = mapping) +
    #geom_point(alpha = 0.3, size = 0.5) +
    #geom_density_2d(color = "orange")
    geom_bin2d() +
    scale_fill_distiller(palette='RdBu', trans='log10') +
    scale_x_continuous(breaks=scales::pretty_breaks())
}
GGally::ggpairs(df, upper = NULL,
                lower = list(continuous = custom_func))

## within genomic compartments

annotsgene <- c("hg19_genes_promoters")
annotations_genes = build_annotations(genome = 'hg19', annotations = annotsgene)
annotations_genes <- unique(annotations_genes) #unique transcripts, not genes
annotations_genes <- annotations_genes[!duplicated(annotations_genes$gene_id)]

annotsgene <- c("hg19_cpg_islands")
annotations_genes = build_annotations(genome = 'hg19', annotations = annotsgene)

grab_sites <- function(regions, gr){
hits <- subjectHits(findOverlaps(regions, gr))
dfsub <- df[hits,]
}

df_proms <- grab_sites(annotations_genes, gr)
df_proms <- df_proms[rowSums(df_proms) > 0,]
e3 <- tidyr::gather(df_proms, Contrast, Change)
e3$Contrast <- factor(e3$Contrast, levels = c("SSA.P", "cADN", "Age", "Segment"))

ggplot(e3) +
  geom_violin(aes(x=Contrast, y=abs(Change)) , fill = "purple") +
  geom_boxplot(aes(x=Contrast, y=abs(Change)), color = "grey", alpha = 0) +
  scale_y_continuous(trans='sqrt') +
  theme_bw() +
  theme(strip.background = element_rect(colour = "black", fill = "white"),
        text = element_text(size = 15),
        legend.position = "none")

GGally::ggpairs(df_proms, upper = NULL,
                lower = list(continuous = custom_func))

