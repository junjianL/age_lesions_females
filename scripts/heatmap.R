## Draw heatmap
# Dec 17 2019

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ComplexHeatmap)
  library(plyranges)
  library(dplyr)
  library(viridis)
  library(annotatr)
})


#### Supervised-ish heatmap, age DMRs ####
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


#### load combined object
load("data/bsseqCombined.RData")

#Get meth table
gr <- rowRanges(bsCombined)
cov <- getCoverage(bsCombined, type = "Cov")
meth <- getCoverage(bsCombined, type = "M")
meth_vals <- meth /cov
colnames(meth_vals) <- paste0(colData(bsCombined)$patient,".",colData(bsCombined)$lesion)
mcols(gr) <- meth_vals
seqlevels(gr) <- paste0("chr",seqlevels(gr))


#summarize meth per region, per sample
hits <- findOverlaps(DMRsage_annot[DMRsage_annot$state == "hyper"], gr)
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
colData(bsCombined)$state <- ifelse(is.na(colData(bsCombined)$state), 
                                        "Normal", colData(bsCombined)$state)

idx <- colData(bsCombined)$state == "Normal" #select all normals
agg <- as.matrix(gr_dmr)[,-1][1:200,idx] #choose number of DMRs
#methsTR <- asin(2*agg-1)

scagg <- scale_exprs(agg)

#Colors
col <- RColorBrewer::brewer.pal(n = 9, name = "YlGnBu")
col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c(col[1], col[5], col[9]))
#col_fun <- viridis(50)
col_anot <- RColorBrewer::brewer.pal(n = 9, name = "Set1")

#change some annotations
seg <- ifelse(colData(bsCombined)$segment == "C","cecum", colData(bsCombined)$segment)
seg <- ifelse(seg == "A","ascend", seg)

#Annotation
column_ha <- HeatmapAnnotation(Age_group = colData(bsCombined)$age_group[idx], 
                              Segment = seg[idx],
                              Tissue = colData(bsCombined)$lesion[idx],
                              col = list(Age_group = c("young" = col_anot[1], 
                                                       "old"= col_anot[2]),
                                         Segment = c("cecum" = col_anot[3], 
                                                     "sigmoid" = col_anot[4],
                                                     "ascend" = col_anot[5]),
                                         Tissue = c("cecum_old" = col_anot[6],
                                                    "cecum_young" = col_anot[7],
                                                    "Normal_Adenoma"= col_anot[1],
                                                    "Normal_SSA"= col_anot[2],
                                                    "sigmoid_old"=col_anot[9],
                                                    "sigmoid_young"=col_anot[8])))

#Plot
pdf("figures/heatmap_agedmrs_plusnormals_mycol.pdf")
Heatmap(scagg, 
        na_col = "white",
        #column_split = metadata$age_group,
        top_annotation = column_ha,
        clustering_distance_columns = "pearson",
        clustering_method_columns = "complete",
        col = col_fun,
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        show_row_dend = FALSE,
        row_title = "Hyper DMRs", 
        column_title = "Samples",
        column_title_side = "bottom",
        heatmap_legend_param = list(title = "beta_value",
                                    grid_height = unit(1, "cm"),
                                    grid_width = unit(0.5, "cm"),
                                    labels_gp = gpar(fontsize = 10)),
        row_names_gp = gpar(fontsize = 12))
dev.off()


#### most variable promoters? ####

annotsgene <- c("hg19_genes_promoters")
annotations_genes = build_annotations(genome = 'hg19', annotations = annotsgene)
annotations_genes <- unique(annotations_genes)


#summarize meth per region, per sample
hits <- findOverlaps(annotations_genes, gr)
gr$promoter <- NA
gr[subjectHits(hits)]$promoter <- queryHits(hits)

gr_sub <- gr[!is.na(gr$promoter)]

#use plyranges
gr_dmr <- gr_sub %>% 
  group_by(promoter) %>% 
  summarise_at(
    colnames(meth_vals), mean, na.rm=TRUE
  )

#Choose most variable promoters with MAD
# med.att <- apply(as.matrix(gr_dmr)[,-1], 1, median)
# mad <- sweep(as.matrix(gr_dmr)[,-1], 1, med.att, FUN="-")
# mads <- apply(mad, 1, function(x) median(abs(x)))

madr <- apply(as.matrix(gr_dmr)[,-1][,idx], 1, mad)
madr <- rowVars(as.matrix(gr_dmr)[,-1][,idx])
o <- order(madr, decreasing = TRUE)
agg <- as.matrix(gr_dmr)[o,-1][1:2000,idx]
methsTR <- asin(2*agg-1)

scagg <- scale_exprs(agg)
  
pdf("figures/heatmap_2000proms_plusnormals_mycol_scaled.pdf")
Heatmap(
        #scagg, 
        methsTR,
        na_col = "white",
        #column_split = metadata$age_group,
        top_annotation = column_ha,
        clustering_distance_columns = "pearson",
        clustering_method_columns = "complete",
        col = col_fun,
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        show_row_dend = FALSE,
        row_title = "Top 200 promoters (MAD)", 
        column_title = "Samples",
        column_title_side = "bottom",
        heatmap_legend_param = list(title = "mean beta value",
                                    grid_height = unit(1, "cm"),
                                    grid_width = unit(0.5, "cm"),
                                    labels_gp = gpar(fontsize = 10)),
        row_names_gp = gpar(fontsize = 12))
dev.off()


scale_exprs <- function(x, margin = 1, q = 0.01) {
  if (!is(x, "matrix")) x <- as.matrix(x)
  qs <- c(rowQuantiles, colQuantiles)[[margin]]
  qs <- qs(x, probs = c(q, 1-q))
  qs <- matrix(qs, ncol = 2)
  x <- switch(margin,
              "1" = (x - qs[, 1]) / (qs[, 2] - qs[, 1]),
              "2" = t((t(x) - qs[, 1]) / (qs[, 2] - qs[, 1])))
  x[x < 0 | is.na(x)] <- 0
  x[x > 1] <- 1
  return(x)
}

