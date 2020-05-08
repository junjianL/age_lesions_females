#!/usr/bin/env Rscript

library(ggplot2)
library(tibble)
library(cowplot)
library(GenomicRanges)
library(annotatr)


load("data/DMRs_lesions_3.RData")
load("data/DMRs_age_final.RData")


#dirty filtering
DMRsles_annot$state <- ifelse(DMRsles_annot$beta > 0 & DMRsles_annot$qval < 0.05 , 
                              "hyper", "none")
DMRsles_annot$state <- ifelse(DMRsles_annot$beta < 0 & DMRsles_annot$qval < 0.05 , 
                              "hypo", DMRsles_annot$state)

DMRsage_annot$state <- ifelse(DMRsage_annot$beta > 0 & DMRsage_annot$qval < 0.05 , 
                              "hyper", "none")
DMRsage_annot$state <- ifelse(DMRsage_annot$beta < 0 & DMRsage_annot$qval < 0.05 , 
                              "hypo", DMRsage_annot$state)

#### annotatr annotation plots ####

annotscpg = c("hg19_cpg_islands", "hg19_cpg_shores",
              "hg19_cpg_shelves")

annotations = build_annotations(genome = 'hg19', annotations = annotscpg)

#genes
annotsgene <- c("hg19_genes_promoters", "hg19_genes_3UTRs", "hg19_genes_introns", 
                "hg19_genes_exons", "hg19_genes_5UTRs", "hg19_genes_cds",
                "hg19_genes_1to5kb")

annotations_genes = build_annotations(genome = 'hg19', annotations = annotsgene)

pdf("figures/annotation_plots_simple.pdf", width = 8, height = 4)

dm_annotated = annotate_regions(
  regions = DMRsles_annot,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)


a <- plot_categorical(
  annotated_regions = dm_annotated, x='state',
  fill='annot.type',
  position='fill',
  plot_title = 'DMRs in lesions by CpG annotation proportions',
  legend_title = 'Annotations',
  x_label = 'Direction',
  y_label = 'Proportion')

a2 <- plot_categorical(
  annotated_regions = dm_annotated, x='state', 
  fill='annot.type',
  position='stack',
  plot_title = 'DMRs in lesions by CpG annotation counts',
  legend_title = 'Annotations',
  x_label = 'Direction',
  y_label = 'Count')

dm_annotated = annotate_regions(
  regions = DMRsage_annot,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

b <- plot_categorical(
  annotated_regions = dm_annotated, x='state',
  fill='annot.type',
  position='fill',
  plot_title = 'DMRs in age by CpG annotation proportions',
  legend_title = 'Annotations',
  x_label = 'Direction',
  y_label = 'Proportion')

b2 <- plot_categorical(
  annotated_regions = dm_annotated, x='state', 
  fill='annot.type',
  position='stack',
  plot_title = 'DMRs in age by CpG annotation counts',
  legend_title = 'Annotations',
  x_label = 'Direction',
  y_label = 'Count')

plot_grid(a2, b2)
plot_grid(a, b)

dev.off()


#### Annotation DMRs ####

# draw_annot.perc <- function(res){
#   inter <- length(res$gene[grep("^intergenic\\.NA$", res$gene, perl=T)]) 
#   proc.inter <- inter /  length(res$gene) * 100
#   intra <- ((length(res$gene) - inter )/ length(res$gene)) * 100
#   d <- data.frame(location = c("intergenic", "intragenic"), DMRs = c(proc.inter, intra))
#   return(d)
# }

#converge info from hypo and hyper
build_plot <- function(tab1, tab2, orientation, orientation2, drawfunc){
  d1 <- cbind(drawfunc(tab1[tab1$state == orientation,]), state = orientation, comparison = "old-young")
  d2 <- cbind(drawfunc(tab2[tab2$state == orientation,]), state = orientation, comparison = "lesion-normal")
  d3 <- cbind(drawfunc(tab1[tab1$state == orientation2,]), state = orientation2, comparison = "old-young")
  d4 <- cbind(drawfunc(tab2[tab2$state == orientation2,]), state = orientation2, comparison = "lesion-normal")
  
  d <- rbind(d1,d2,d3,d4)
  
  p1 <- ggplot(data=d, aes(x=location, y=number.DMRs, fill = comparison)) + #change number for percentage
    geom_bar(stat="identity", position=position_dodge()) + 
    theme_bw() +
    theme(legend.position="bottom", text = element_text(size=15)) +
    labs(x = "Location", y = "Number of DMRs") +
    facet_grid(~state)
  return(p1)
}

#genic annotation
draw_annot_genic <- function(res){
  inter <- length(res$gene[grep("^intergenic\\.NA$", res$gene, perl=T)]) 
  #proc.inter <- inter /  length(res$gene) * 100
  intra <- length(res$gene) - inter 
  d <- data.frame(location = c("intergenic", "intragenic"), number.DMRs = c(inter, intra))
  return(d)
}

#intragenic annotation
draw_annot_intragenic <- function(res){
  UTR3 <- grep("3UTR", res$gene)
  res_no3UTR <- res[-UTR3,]
  proms <- grep("(promoter|1to5kb|5UTR)", res_no3UTR$gene, perl =T)
  res_no3UTR_noproms <- res_no3UTR[-proms,]
  gene_body <- grep("(exon|intron|CDS)", res_no3UTR_noproms$gene, perl =T)
  #res_no3UTR_nogeneb_noproms <- res_no3UTR_noproms[-gene_body,]
  d <- tibble(location = factor(c("Regulatory region", "Gene body", "3UTR")),
              number.DMRs = c(length(proms), length(gene_body), length(UTR3)))
  d$location <- factor(d$location, c("Regulatory region", "Gene body", "3UTR"))
  return(d)
}

#CpG related annotation
draw_annot_cpg_general <- function(res){
  inter <- length(res$CpG[grep("^inter:[0-9]+$", res$CpG, perl = T)])
  intra <- length(res$CpG) - inter
  d <- data.frame(location = c("Outside SSI", "SSI"), number.DMRs = c(inter, intra))
  return(d)
}

#CpG island related annotation
draw_annot_cpg <- function(res){
  islands <- grep("island", res$CpG)
  res_noisl <- res[-islands,]
  shore <- grep("shore", res_noisl$CpG)
  res_noisl_noshore <- res_noisl[-shore,]
  shelf <- grep("shelf", res_noisl_noshore$CpG)
  res_noisl_noshore_noshelf <- res_noisl_noshore[-shelf,]
  d <- tibble(location = factor(c("CpG island", "CpG shore", "CpG shelf")),
              number.DMRs = c(length(islands), length(shore), length(shelf)))
  d$location <- factor(d$location, c("CpG island", "CpG shore", "CpG shelf"))
  return(d)
  
}

pdf("figures/annotation_plots.pdf", width = 8, height = 5)
build_plot(DMRsage_annot, DMRsles_annot, "hyper", "hypo", draw_annot_genic)
build_plot(DMRsage_annot, DMRsles_annot, "hyper", "hypo", draw_annot_intragenic)
build_plot(DMRsage_annot, DMRsles_annot, "hyper", "hypo", draw_annot_cpg_general)
build_plot(DMRsage_annot, DMRsles_annot, "hyper", "hypo", draw_annot_cpg)

dev.off()
