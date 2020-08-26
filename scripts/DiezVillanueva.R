## Public data
# August 4 2020

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(stringr)
  library(RColorBrewer)
  library(minfi)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(ComplexHeatmap)
  library(dplyr)
  library(plyranges)
})

## DiezVillanueva data

#### re-format metadata ####
metadatafull <- readr::read_delim("data/public_data/DiezVillanueva_2020/metadata.txt",
                              delim = "\t",
                              col_names = FALSE)
metadatafull <- t(metadatafull)[-1,-1]
metadatafull <- as.data.frame(metadatafull)


#there's a column missing in the healthy samples, so easiest is to split early
#and then recombine

#split healthy
metadata <- metadatafull[,-6]

idxh <- grepl("_M", metadata$V1)
#hsamps <- metadata$V1[idxh]
metadata <- metadata[idxh,]

colnames(metadata) <- c("sample", "tissue", "gender", "age", "location")

#split crc
metadatacrc <- metadatafull[!idxh,]
head(metadatacrc)
metadatacrc <- metadatacrc[,-3]
colnames(metadatacrc) <- c("sample", "tissue", "gender", "age", "location")

#join
metadata <- rbind(metadata, metadatacrc)

#clean var types
head(metadata)
metadata$age <- as.integer(gsub("age: ","",metadata$age))
metadata$age_group <- "mid-age"
metadata$age_group[metadata$age > 70] <- "old"
metadata$age_group[metadata$age <= 40] <- "young"
metadata$age_group <- factor(metadata$age_group, levels = c("young", "mid-age", "old"))
table(metadata$age_group)
#only 1 sample under 40

metadata$tissue <- factor(gsub("tissue: ","",metadata$tissue), levels=c("Mucosa", "Normal", "Tumor"))
metadata$location <- factor(gsub("location: ","",metadata$location), levels = c("Left", "Right"))

#### Get GEO data ####
# mset <- getGenomicRatioSetFromGEO("GSE131013")
# Error in getGenomicRatioSetFromGEO("GSE131013") : 
#   No rowname matches. 'rownames' need to match IlluminaHumanMethylation450k probe names.
# doesnt have the beta values in the SOFT files, was loaded separately

# Import normalized matrix 

betas <- read.table("data/public_data/DiezVillanueva_2020/GSE131013_normalized_matrix.txt.gz", 
                    header = TRUE, row.names = 1)

#separate tables
pvals <- betas[,grepl("Detection.Pval", colnames(betas))]
colnames(pvals) <- gsub(".Detection.Pval","",colnames(pvals))

betas <- betas[,!grepl("Detection.Pval", colnames(betas))]

#reorder like metadata 
pvals <- pvals[,metadata$sample]
betas <- betas[,metadata$sample]

# get the 450k annotation data
ann450k <- getLocations(IlluminaHumanMethylation450kanno.ilmn12.hg19,
                        orderByLocation = TRUE)

ind2 <- match(rownames(betas),names(ann450k))

# remove non-matching probes
rem <- which(is.na(ind2))
betas <- betas[-rem,] 
pvals <- pvals[-rem,] 
ind2 <- ind2[-rem] 

dim(betas) == dim(pvals)

#build GRset
diez <- GenomicRatioSet(ann450k[ind2], Beta = betas, M = NULL,
                        CN = NULL,
                        annotation = c(array = "IlluminaHumanMethylation450k",
                                       annotation = "ilmn12.hg19"),
                        colData = metadata)
colnames(diez) <- colnames(betas)
assay(diez, "pvalue") <- pvals

diez$dataset <- "DiezVillanueva"
saveRDS(diez, "data/public_data/diez.rds")

# choose females
diezfem <- diez[,diez$gender == "gender: Female"]

# remove low detection pvalues across all samples
idx <- rowSums(assay(diezfem, "pvalue") < 0.05) == ncol(diezfem)
diezfem <- diezfem[idx,]

#do some plots
pal <- brewer.pal(8,"Dark2")

#by tissue
limma::plotMDS(getM(diezfem), top=10000, gene.selection="common",
               col=pal[factor(diezfem$tissue)])
legend("top", legend=levels(factor(diezfem$tissue)), text.col=pal,
       bg="white", cex=0.7)


#### draw heatmap with selected markers ####

load("data/rdata/unique_lesions_filt.RData")

# get probes in selected markers
meth_vals <- as.matrix(getBeta(diezfem))

draw_hm <- function(obj, regions, numregs){
  gr <- rowRanges(obj)
  mcols(gr) <- meth_vals
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
  madr <- rowVars(as.matrix(gr_dmr)[,-1])
  o <- order(madr, decreasing = TRUE)
  score <- as.matrix(gr_dmr)[o,-1][1:numregs,]
  
  return(score)
  
  #Colors
  col <- RColorBrewer::brewer.pal(n = 9, name = "YlGnBu")
  col_anot <- RColorBrewer::brewer.pal(n = 9, name = "Set1")
  purples <- RColorBrewer::brewer.pal(n = 3, name = "Purples")
  greens <- RColorBrewer::brewer.pal(n = 3, name = "Greens")
  pinks <- RColorBrewer::brewer.pal(n = 4, name = "RdPu")[c(1,3:4)]
  
  #hm colors
  col_fun <- circlize::colorRamp2(c(0,0.2,1), c(col[9], col[7], col_anot[6]))
  
  #annot colors
  col_age <- purples[1:nlevels(obj$age_group)]
  names(col_age) <- levels(obj$age_group)
  
  col_tis <- greens[1:nlevels(obj$tissue)]
  names(col_tis) <- levels(obj$tissue)
  
  col_seg <- pinks[1:nlevels(obj$location)]
  names(col_seg) <- levels(obj$location)
  
  #column annot
  column_ha <- HeatmapAnnotation(Age = obj$age_group, 
                                 Tissue = obj$tissue,
                                 Location = obj$location,
                                 col = list(Age = col_age,
                                            Tissue = col_tis,
                                            Location = col_seg), 
                                 gp = gpar(col = "black"))
  
  #remove names
  rownames(score) <- colnames(score) <- NULL
  
  #Plot
  hm <- Heatmap(score, 
                use_raster = TRUE,
                na_col = "white",
                column_split = obj$tissue,
                top_annotation = column_ha,
                col = col_fun,
                #row_km = 2, 
                clustering_distance_columns = "spearman",
                cluster_columns = TRUE,
                show_row_dend = FALSE,
                show_column_dend = TRUE,
                cluster_column_slices = FALSE,
                #left_annotation = row_ha, ## remove with 3 groups
                row_title = "selected markers", 
                column_title = "Samples",
                column_title_side = "bottom",
                column_names_gp = gpar(fontsize = 8),
                heatmap_legend_param = list(title = "Beta",
                                            title_position = "lefttop-rot",
                                            grid_height = unit(1, "cm"),
                                            grid_width = unit(0.5, "cm")))
  #return(hm)
}

draw_hm(diezfem, sub_unique, 2000)

## non-selected regions
library(annotatr)
annotsgene <- c("hg19_cpg_islands")
annotations_genes <- build_annotations(genome = 'hg19', annotations = annotsgene)


#### combine all markers into single signature ####

combine_marks <- function(obj, regions){
  gr <- rowRanges(obj)
  mcols(gr) <- meth_vals
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
  score <- as.matrix(gr_dmr)[,-1]
  mean_betas <- colMedians(score, na.rm = TRUE)
}
mean_betas <- combine_marks(diezfem, sub_uniqueannot)

#set truth and draw ROC
truth <- ifelse(grepl("_T",colnames(meth_vals)), 1,0)
pROC::roc(truth, mean_betas, direction = "<", plot = TRUE,
          main= "Combined selected markers, Diez-Villanueva 2020",
          percent=TRUE, print.auc=TRUE, print.thres = "best")

