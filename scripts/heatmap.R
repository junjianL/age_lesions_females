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

#### Function to get matrix, to get annotations, calculate seriation, draw HM ####

hm_build <- function(regions, sampleidx, numregs=2000, mostvar = FALSE, 
                     split = "age_group", includeseg = FALSE,
                     title = "Most variable 2000 CpG promoters", numgroups = 1){
  
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
  agg[is.na(agg)] <- 0
  agg[is.nan(agg)] <- 0
  agg[is.infinite(agg)] <- 0
  
  #row annotation (meth change)
  if (numgroups == 2) {
  condit <- colData(bsCombined)[,split][idx]
  conds <- as.character(unique(condit))
  meth_change <- rowMeans(agg[,condit == conds[2]]) - rowMeans(agg[,condit == conds[1]])
  }
  
  # ord <- seriation::get_order(seriation::seriate(dist(t(agg)), 
  #                                                method = "MDS_angle"))
  #                                                #margin = 2))
  #agg <- agg[,ord]
  
  #Colors
  col <- RColorBrewer::brewer.pal(n = 9, name = "YlGnBu")
  col_anot <- RColorBrewer::brewer.pal(n = 9, name = "Set1")
  purples <- RColorBrewer::brewer.pal(n = 3, name = "Purples")
  oranges <- RColorBrewer::brewer.pal(n = 3, name = "Oranges")
  greens <- RColorBrewer::brewer.pal(n = 3, name = "Greens")
  pinks <- RColorBrewer::brewer.pal(n = 4, name = "RdPu")[c(1,3:4)]
  
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

  column_ha <- HeatmapAnnotation(Age = colData(bsCombined)$age_group[idx], #[ord] 
                                 Diagnosis = colData(bsCombined)$diagnosis[idx], #[ord]
                                 Tissue = colData(bsCombined)$tissue[idx], #[ord]
                                 col = list(Age = col_age,
                                            Diagnosis = col_diag,
                                            Tissue = col_tis
                                 ), 
                                 gp = gpar(col = "black"))
  
  if (numgroups == 2) {
  row_ha <- rowAnnotation("Change" = anno_lines(meth_change, 
                                              smooth = FALSE, #loess
                                              add_points = FALSE,
                                              axis_param = list(direction = "reverse",
                                                                gp = gpar(fontsize = 5)
                                                                )
                                              ))
  }

  if(includeseg) {
    column_ha <- HeatmapAnnotation(Age = colData(bsCombined)$age_group[idx],#[ord]
                                   Segment = colData(bsCombined)$seg[idx], #[ord]
                                   Diagnosis = colData(bsCombined)$diagnosis[idx], #[ord]
                                   Tissue = colData(bsCombined)$tissue[idx], #[ord]
                                   col = list(Age = col_age,
                                              Diagnosis = col_diag,
                                              Tissue = col_tis,
                                              Segment = col_seg
                                   ), gp = gpar(col = "black"))
  }
  
  #Plot
  hm <- Heatmap(agg, 
          na_col = "white",
          column_split = colData(bsCombined)[,split][idx], #[ord]
          top_annotation = column_ha,
          col = col_fun,
          row_km = 2, 
          clustering_distance_columns = "spearman",
          cluster_columns = TRUE,
          show_row_dend = FALSE,
          show_column_dend = TRUE,
          cluster_column_slices = FALSE,
          #left_annotation = row_ha, ## remove with 3 groups
          row_title = title, 
          column_title = "Samples",
          column_title_side = "bottom",
          column_names_gp = gpar(fontsize = 8),
          heatmap_legend_param = list(title = "mean beta value",
                                      title_position = "lefttop-rot",
                                      grid_height = unit(1, "cm"),
                                      grid_width = unit(0.5, "cm")))
                                      #labels_gp = gpar(fontsize = 6)))
  
  return(hm)
}
  
#### Plotting supervised HM ####

#Get top regions
load("data/DMRs_age_final.RData")

#filtering
DMRsage_annot$state <- ifelse(DMRsage_annot$beta > 0 & DMRsage_annot$qval < 0.05 , 
                              "hyper", "none")
DMRsage_annot$state <- ifelse(DMRsage_annot$beta < 0 & DMRsage_annot$qval < 0.05 , 
                              "hypo", DMRsage_annot$state)

idx <- colData(bsCombined)$state == "Normal"

hm_build(DMRsage_annot[DMRsage_annot$state == "hyper"], idx, numregs= 1000,
         title = "Top 1000 age-DMRs", split = "age_group", includeseg = TRUE)
dev.off()
#### most variable promoters? ####

annotsgene <- c("hg19_genes_promoters")

annotsgene <- c("hg19_cpg_islands")

annotations_genes = build_annotations(genome = 'hg19', annotations = annotsgene)
annotations_genes <- unique(annotations_genes) #unique transcripts, not genes
annotations_genes <- annotations_genes[!duplicated(annotations_genes$gene_id)]

idx <- colData(bsCombined)$tissue == "normal mucosa" #select all normals
idx <- colData(bsCombined)$state == "Normal" #selet all healthy fems
idx <- rep(TRUE, ncol(bsCombined)) # all samples

hm_build(annotations_genes, idx, numregs= 1000, split = "tissue", includeseg = TRUE,
         title = "1000 most variable promoters", mostvar = TRUE)


hm_build(annotations_genes, idx, numregs= 1000, split = "age_group", numgroups = 1,
         title = "1000 most variable promoters", mostvar = TRUE, includeseg = TRUE)

#### horvath age probes ####
annot450k <- readr::read_csv("data/HumanMethylation450_15017482_v1-2_edited.csv")

#probes from hovarth
probes <- read.csv("data/13059_2013_3156_MOESM3_ESM.csv", header = TRUE)

#get buil 37 locations from annot450k
cidx <- match(probes$CpGmarker, annot450k$Name)
cidx <- cidx[!is.na(cidx)]
chr <- annot450k$CHR[idx]
pos <- annot450k$MAPINFO[idx]

idx <- colData(bsCombined)$tissue == "healthy" #select all normals
idx <- colData(bsCombined)$state == "Normal"

probesgr <- GRanges(paste0("chr",chr), IRanges(start = pos-1, end = pos))

hm_build(probesgr, idx, numregs= 317, split = "age_group",
         title = "aging clock sites", includeseg = TRUE)

