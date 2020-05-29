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


seg <- ifelse(colData(bsCombined)$segment == "C","cecum", colData(bsCombined)$segment)
seg <- ifelse(seg == "A","ascend", seg)

age <- ifelse(colData(bsCombined)$age < 40, "<40", "41-70")
age <- ifelse(colData(bsCombined)$age > 70, ">70", age)
colData(bsCombined)$age_group <- age

#summarize meth per region, per sample <<<<<<-------------------
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

#### Get matrix for all normals ####
idx <- colData(bsCombined)$state == "Normal" #select all normals
#or
idx <- rep(TRUE, ncol(bsCombined))
agg <- as.matrix(gr_dmr)[,-1][1:200,idx] #choose number of DMRs

#methsTR <- asin(2*agg-1)
#scagg <- scale_exprs(agg)

#### Plotting supervised HM ####
#Colors
#my prefered color scheme
col <- RColorBrewer::brewer.pal(n = 9, name = "YlGnBu")
#col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c(col[1], col[5], col[9]))

#annotation colors
col_anot <- RColorBrewer::brewer.pal(n = 9, name = "Set1")
col_anot2 <- RColorBrewer::brewer.pal(n = 3, name = "Set2")

#giancarlo schemes
#col_fun <- viridis(50) #1
col_fun <- circlize::colorRamp2(c(0,0.2,1), c(col[9], col[7], col_anot[6]))


#HM Annotation
column_ha <- HeatmapAnnotation(Age_group = colData(bsCombined)$age_group[idx], 
                              Segment = seg[idx],
                              Condition = colData(bsCombined)$state[idx],
                              col = list(Age_group = c("<40" = col_anot2[1], 
                                                       "41-70"= col_anot2[2],
                                                       ">70"= col_anot2[3]),
                                         Segment = c("cecum" = col_anot[3], 
                                                     "sigmoid" = col_anot[4],
                                                     "ascend" = col_anot[5]),
                                         # Tissue = c("cecum_old" = col_anot[6],
                                         #            "cecum_young" = col_anot[7],
                                         #            "Normal_Adenoma"= col_anot[1],
                                         #            "Normal_SSA"= col_anot[2],
                                         #            "sigmoid_old"=col_anot[9],
                                         #            "sigmoid_young"=col_anot[8]))
                                         Condition = c("Normal" = col_anot[6],
                                                       "Normal_SSA" = col_anot[7],
                                                       "Normal_cADN" = col_anot[1],
                                                       "SSA" = col_anot[9],
                                                       "cADN" = col_anot[8])
                                         ),
                              gp = gpar(col = "black"))

# only normal annotation
column_ha <- HeatmapAnnotation(Age_group = colData(bsCombined)$age_group[idx], 
                               Segment = seg[idx],
                               Condition = colData(bsCombined)$state[idx],
                               col = list(Age_group = c("<40" = col_anot2[1], 
                                                        "41-70"= col_anot2[2],
                                                        ">70"= col_anot2[3]),
                                          Segment = c("cecum" = col_anot[3], 
                                                      "sigmoid" = col_anot[4],
                                                      "ascend" = col_anot[5]),
                                          # Tissue = c("cecum_old" = col_anot[6],
                                          #            "cecum_young" = col_anot[7],
                                          #            "Normal_Adenoma"= col_anot[1],
                                          #            "Normal_SSA"= col_anot[2],
                                          #            "sigmoid_old"=col_anot[9],
                                          #            "sigmoid_young"=col_anot[8]))
                                          Condition = c("Normal" = col_anot[6],
                                                        "Normal_SSA" = col_anot[7],
                                                        "Normal_cADN" = col_anot[1]
                                                        #"SSA" = col_anot[9],
                                                        #"cADN" = col_anot[8]
                                                        )
                               ),
                               gp = gpar(col = "black"))

#Plot
pdf("figures/heatmap_agedmrs_allsamples.pdf")
Heatmap(agg, 
        na_col = "white",
        #column_split = colData(bsCombined)$age_group[idx],
        #column_split = seg[idx],
        top_annotation = column_ha,
        clustering_distance_columns = "spearman",
        clustering_method_columns = "complete",
        col = col_fun,
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        show_row_dend = FALSE,
        row_title = "age related hyper-DMRs (200)", 
        column_title = "Samples",
        column_title_side = "bottom",
        column_names_gp = gpar(fontsize = 8),
        heatmap_legend_param = list(title = "mean beta",
                                    grid_height = unit(1, "cm"),
                                    grid_width = unit(0.5, "cm"),
                                    labels_gp = gpar(fontsize = 10)),
        row_names_gp = gpar(fontsize = 12))
dev.off()

##### Only use age samples ####
idx <- colData(bsCombined)$segment %in% c("cecum","sigmoid") #select all normals
agg <- as.matrix(gr_dmr)[,-1][1:200,idx] #choose number of DMRs

column_ha <- HeatmapAnnotation(Age_group = colData(bsCombined)$age_group[idx], 
                               Segment = seg[idx],
                               col = list(Age_group = c("<40" = col_anot2[1], 
                                                        "41-70"= col_anot2[2],
                                                        ">70"= col_anot2[3]),
                                          Segment = c("cecum" = col_anot[3], 
                                                      "sigmoid" = col_anot[4],
                                                      "ascend" = col_anot[5])
                                          ))


pdf("figures/heatmap_agedmrs_splitage.pdf")
Heatmap(agg, 
        na_col = "white",
        column_split = colData(bsCombined)$age_group[idx],
        #column_split = seg[idx],
        top_annotation = column_ha,
        clustering_distance_columns = "spearman",
        clustering_method_columns = "complete",
        col = col_fun,
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        show_row_dend = FALSE,
        row_title = "age related hyper-DMRs (top 200)", 
        column_title = "Samples",
        column_title_side = "bottom",
        heatmap_legend_param = list(title = "beta",
                                    grid_height = unit(1, "cm"),
                                    grid_width = unit(0.5, "cm"),
                                    labels_gp = gpar(fontsize = 10)),
        row_names_gp = gpar(fontsize = 12))
dev.off()

#### most variable promoters? ####

annotsgene <- c("hg19_genes_promoters")
#annotsgene <- c("hg19_cpg_islands")
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

idx <- grepl("Normal", colData(bsCombined)$state)
idx <- rep(TRUE, ncol(bsCombined))

#madr <- apply(as.matrix(gr_dmr)[,-1][,idx], 1, mad)
madr <- rowVars(as.matrix(gr_dmr)[,-1][,idx])
o <- order(madr, decreasing = TRUE)
agg <- as.matrix(gr_dmr)[o,-1][1:2000,idx]
#methsTR <- asin(2*agg-1)

pdf("figures/heatmap_2000proms_allsamps_splitcond.pdf")
Heatmap(agg, 
        #methsTR,
        na_col = "white",
        #column_split = colData(bsCombined)$age_group[idx],
        #column_split = seg[idx],
        column_split = colData(bsCombined)$state[idx],
        top_annotation = column_ha,
        clustering_distance_columns = "spearman",
        clustering_method_columns = "complete",
        col = col_fun,
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        row_title = "Most variable 2000 CpG promoters", 
        column_title = "Samples",
        column_title_side = "bottom",
        column_names_gp = gpar(fontsize = 8),
        heatmap_legend_param = list(title = "mean beta value",
                                    grid_height = unit(1, "cm"),
                                    grid_width = unit(0.5, "cm"),
                                    labels_gp = gpar(fontsize = 10)),
        row_names_gp = gpar(fontsize = 12))
dev.off()


## only age samples
idx <- colData(bsCombined)$segment %in% c("sigmoid", "cecum")

#madr <- apply(as.matrix(gr_dmr)[,-1][,idx], 1, mad)
madr <- rowVars(as.matrix(gr_dmr)[,-1][,idx])
o <- order(madr, decreasing = TRUE)
agg <- as.matrix(gr_dmr)[o,-1][1:200,idx]
methsTR <- asin(2*agg-1)

column_ha <- HeatmapAnnotation(Age_group = colData(bsCombined)$age_group[idx], 
                               Segment = seg[idx],
                               col = list(Age_group = c("<40" = col_anot2[1], 
                                                        "41-70"= col_anot2[2],
                                                        ">70"= col_anot2[3]),
                                          Segment = c("cecum" = col_anot[3], 
                                                      "sigmoid" = col_anot[4],
                                                      "ascend" = col_anot[5])
                               ))

pdf("figures/heatmap_200cpgisles_splitsegment.pdf")
Heatmap(
  agg, 
  #methsTR,
  na_col = "white",
  #column_split = colData(bsCombined)$age_group[idx],
  column_split = seg[idx],
  top_annotation = column_ha,
  clustering_distance_columns = "spearman",
  clustering_method_columns = "complete",
  col = col_fun,
  cluster_rows = TRUE, 
  cluster_columns = TRUE,
  show_row_dend = FALSE,
  row_title = "Most variable 200 CpG Islands", 
  column_title = "Samples",
  column_title_side = "bottom",
  heatmap_legend_param = list(title = "mean beta value",
                              grid_height = unit(1, "cm"),
                              grid_width = unit(0.5, "cm"),
                              labels_gp = gpar(fontsize = 10)),
  row_names_gp = gpar(fontsize = 12))
dev.off()
