## Draw heatmap
# Dec 17 2019

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ComplexHeatmap)
  library(plyranges)
  library(dplyr)
  library(viridis)
  library(annotatr)
  library(bsseq)
})


#### load combined object ####
load("data/bsseqCombined.RData")

#Get meth table
gr <- rowRanges(bsCombined)
cov <- getCoverage(bsCombined, type = "Cov")
meth <- getCoverage(bsCombined, type = "M")
meth_vals <- meth /cov
colnames(meth_vals) <- paste0(colData(bsCombined)$patient,".",colData(bsCombined)$lesion)
mcols(gr) <- meth_vals
seqlevels(gr) <- paste0("chr",seqlevels(gr))

#Fix some annotations
colData(bsCombined)$state <- ifelse(is.na(colData(bsCombined)$state), 
                                    "Normal", colData(bsCombined)$state)

colData(bsCombined)$state <- ifelse(grepl("SSA|Adenoma", colData(bsCombined)$lesion), 
                                    colData(bsCombined)$lesion, colData(bsCombined)$state)

colData(bsCombined)$state <- ifelse(grepl("Adenoma", colData(bsCombined)$state), 
                                    gsub("Adenoma","cADN", colData(bsCombined)$state), 
                                    colData(bsCombined)$state)

diagnosis <- ifelse(grepl("SSA|cADN", colData(bsCombined)$state), 
                                        colData(bsCombined)$state, "healthy")

colData(bsCombined)$diagnosis <- factor(gsub("Normal_","",diagnosis), 
                                        levels = c("healthy", "cADN", "SSA"))

colData(bsCombined)$tissue <- factor(ifelse(grepl("Normal", colData(bsCombined)$state), 
                                     "healthy","pre-lesion"))



seg <- ifelse(colData(bsCombined)$segment == "C","cecum", colData(bsCombined)$segment)
colData(bsCombined)$seg <- factor(ifelse(seg == "A","ascend", seg))


age <- ifelse(colData(bsCombined)$age < 40, "<40", "41-70")
age <- ifelse(colData(bsCombined)$age > 70, ">70", age)
colData(bsCombined)$age_group <- factor(age, levels = c("<40", "41-70", ">70"))

#### Function to get matrix, to get annotations, calculate seriation, draw HM ####

hm_build <- function(regions, sampleidx, numregs=2000, mostvar = FALSE, 
                     annot = "all", split = "age_group", includeseg = FALSE,
                     title = "Most variable 2000 CpG promoters"){
  
  hits <- findOverlaps(regions, gr)
  gr$DMR <- NA
  gr[subjectHits(hits)]$DMR <- queryHits(hits)
  
  gr_sub <- gr[!is.na(gr$DMR)]
  
  #use plyranges
  gr_dmr <- gr_sub %>% 
    group_by(DMR) %>% 
    summarise_at(
      colnames(meth_vals), mean, na.rm=TRUE
    )
  
  #Get matrix
  if(mostvar) {
  madr <- rowVars(as.matrix(gr_dmr)[,-1][,sampleidx])
  o <- order(madr, decreasing = TRUE)
  agg <- as.matrix(gr_dmr)[o,-1][1:numregs,sampleidx]
  } else agg <- as.matrix(gr_dmr)[,-1][1:numregs,sampleidx] 

  #seriation to order samples
  
  ord <- seriation::get_order(seriation::seriate(dist(t(agg)), 
                                                 method = "MDS_angle"))
                                                 #margin = 2))
  agg <- agg[,ord]
  
  #Colors
  col <- RColorBrewer::brewer.pal(n = 9, name = "YlGnBu")
  col_anot <- RColorBrewer::brewer.pal(n = 9, name = "Set1")
  purples <- RColorBrewer::brewer.pal(n = 3, name = "Purples")
  oranges <- RColorBrewer::brewer.pal(n = 3, name = "Oranges")
  greens <- RColorBrewer::brewer.pal(n = 3, name = "Greens")
  pinks <- RColorBrewer::brewer.pal(n = 3, name = "RdPu")
  
  #hm colors
  col_fun <- circlize::colorRamp2(c(0,0.2,1), c(col[9], col[7], col_anot[6]))
  
  #annot colors
  col_age <- purples[1:nlevels(colData(bsCombined)$age_group[idx])]
  names(col_age) <- levels(colData(bsCombined)$age_group[idx])
  
  col_diag <- oranges[1:nlevels(colData(bsCombined)$diagnosis[idx])]
  names(col_diag) <- levels(colData(bsCombined)$diagnosis[idx])
  
  col_tis <- greens[1:nlevels(colData(bsCombined)$tissue[idx])]
  names(col_tis) <- levels(colData(bsCombined)$tissue[idx])
  
  col_seg <- pinks[1:nlevels(colData(bsCombined)$seg[idx])]
  names(col_seg) <- levels(colData(bsCombined)$seg[idx])
  #HM Annotation
  # } 
  # if (annot == "healthy") {
  #   condition <- c("Normal" = col_anot[6])
  # }

  column_ha <- HeatmapAnnotation("Age group" = colData(bsCombined)$age_group[idx][ord], 
                                 Diagnosis = colData(bsCombined)$diagnosis[idx][ord],
                                 Tissue = colData(bsCombined)$tissue[idx][ord],
                                 col = list("Age group" = col_age,
                                            Diagnosis = col_diag,
                                            Tissue = col_tis
                                 ), gp = gpar(col = "black"))
  
  if(includeseg) {
    column_ha <- HeatmapAnnotation("Age group" = colData(bsCombined)$age_group[idx][ord], 
                                   Segment = colData(bsCombined)$seg[idx][ord],
                                   Diagnosis = colData(bsCombined)$diagnosis[idx][ord],
                                   Tissue = colData(bsCombined)$tissue[idx][ord],
                                   col = list("Age group" = col_age,
                                              Diagnosis = col_diag,
                                              Tissue = col_tis,
                                              Segment = col_seg
                                   ), gp = gpar(col = "black"))
  }
  
  #Plot
  hm <- Heatmap(agg, 
          na_col = "white",
          column_split = colData(bsCombined)[,split][idx][ord],
          top_annotation = column_ha,
          col = col_fun,
          row_km = 2, 
          clustering_distance_rows = "spearman",
          cluster_columns = FALSE,
          show_row_dend = FALSE,
          show_column_dend = FALSE,
          row_title = title, 
          column_title = "Samples",
          column_title_side = "bottom",
          column_names_gp = gpar(fontsize = 8),
          heatmap_legend_param = list(title = "mean beta-value",
                                      grid_height = unit(1, "cm"),
                                      grid_width = unit(0.5, "cm"),
                                      labels_gp = gpar(fontsize = 10)),
          row_names_gp = gpar(fontsize = 12))
  
  return(hm)
}
  
#### Plotting supervised HM ####

#Get top regions
load("data/DMRs_lesions_3.RData")
load("data/DMRs_age_final.RData")


#filtering
DMRsles_annot$state <- ifelse(DMRsles_annot$beta > 0 & DMRsles_annot$qval < 0.05 , 
                              "hyper", "none")
DMRsles_annot$state <- ifelse(DMRsles_annot$beta < 0 & DMRsles_annot$qval < 0.05 , 
                              "hypo", DMRsles_annot$state)

DMRsage_annot$state <- ifelse(DMRsage_annot$beta > 0 & DMRsage_annot$qval < 0.05 , 
                              "hyper", "none")
DMRsage_annot$state <- ifelse(DMRsage_annot$beta < 0 & DMRsage_annot$qval < 0.05 , 
                              "hypo", DMRsage_annot$state)

idx <- colData(bsCombined)$state == "Normal" #select all normals
#or
idx <- rep(TRUE, ncol(bsCombined))

hm_build(DMRsage_annot[DMRsage_annot$state == "hyper"], idx, numregs= 100,
         title = "Top 200 age-DMRs")

#### most variable promoters? ####

annotsgene <- c("hg19_genes_promoters")
#annotsgene <- c("hg19_cpg_islands")
annotations_genes = build_annotations(genome = 'hg19', annotations = annotsgene)
annotations_genes <- unique(annotations_genes)

idx <- colData(bsCombined)$tissue == "healthy" #select all normals
idx <- colData(bsCombined)$state == "Normal" #selet all healthy fems
idx <- rep(TRUE, ncol(bsCombined))

hm_build(annotations_genes, idx, numregs= 1000, split = "tissue",
         title = "1000 most variable promoters", mostvar = TRUE)


hm_build(annotations_genes, idx, numregs= 1000, split = "age_group",
         title = "1000 most variable promoters", mostvar = TRUE, includeseg = TRUE)

