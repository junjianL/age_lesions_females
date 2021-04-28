##########################################################
## Code to download data from dataset GSE48684, 
# and draw figure 4A
# August 25 2020
##########################################################

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

## Luo, Grady data
# useful for adenoma samples

#import directly from GEO (this is basically an SE object)
mset <- getGenomicRatioSetFromGEO("GSE48684")

# clean vars
# no age information in GEO

#Stage: characteristics_ch1.4 

mset$location <- mset$`colon region:ch1`
mset$location <- gsub("left", "Left", mset$location)
mset$location <- gsub("right", "Right", mset$location) 
mset$location <- gsub("Retum", "Rectum", mset$location) 
mset$location <- gsub("colon", "Unknown", mset$location) 
#mset$location <- gsub("Proximal|Transverse", "Right", mset$location)
#mset$location <- gsub("Distal", "Left", mset$location)
mset$location <- factor(mset$location, levels = c("Unknown","Right","Proximal",
                                                  "Transverse", "Left","Distal", 
                                                  "Rectum"))

mset$tissue <- factor(mset$description, levels = c("normal-H", "normal-C",
                                                   "adenoma", "cancer"))
mset$dataset <- "Luo"
saveRDS(mset, "data/public_data/luo.rds")
#mset <- readRDS("data/public_data/luo.rds")

#select females
msetfem <- mset[,mset$`gender:ch1` %in% c("female", "Female")]


#### Code not used in paper ####
#draw MDS plots
pal <- brewer.pal(8,"Dark2")

#by tissue
limma::plotMDS(getM(msetfem), top=10000, gene.selection="common",
               col=pal[factor(msetfem$tissue)])
legend("top", legend=levels(factor(msetfem$tissue)), text.col=pal,
       bg="white", cex=0.7)

#by location
limma::plotMDS(getM(msetfem), top=10000, gene.selection="common",
               col=pal[factor(msetfem$location)])
legend("top", legend=levels(factor(msetfem$location)), text.col=pal,
       bg="white", cex=0.7)

#### Draw heatmap ####
meth_vals <- as.matrix(getBeta(msetfem))

draw_hm <- function(obj, regions){
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
  #madr <- rowVars(as.matrix(gr_dmr)[,-1])
  #o <- order(madr, decreasing = TRUE)
  score <- as.matrix(gr_dmr)[,-1]
  print(dim(score))
  
  #return(score)
  
  #Colors
  col <- RColorBrewer::brewer.pal(n = 9, name = "YlGnBu")
  col_anot <- RColorBrewer::brewer.pal(n = 9, name = "Set1")
  #purples <- RColorBrewer::brewer.pal(n = 3, name = "Purples")
  greens <- RColorBrewer::brewer.pal(n = 4, name = "Greens")
  pinks <- RColorBrewer::brewer.pal(n = 7, name = "RdPu")
  
  #hm colors
  col_fun <- circlize::colorRamp2(c(0,0.2,1), c(col[9], col[7], col_anot[6]))
  
  #annot colors
  col_tis <- greens[1:nlevels(obj$tissue)]
  names(col_tis) <- levels(obj$tissue)
  
  col_seg <- pinks[1:nlevels(obj$location)]
  names(col_seg) <- levels(obj$location)
  
  #column annot
  column_ha <- HeatmapAnnotation(#Age = obj$age_group, 
                                 Tissue = obj$tissue,
                                 Location = obj$location,
                                 col = list(#Age = col_age,
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
                clustering_distance_columns = "spearman",
                cluster_columns = TRUE,
                show_row_dend = FALSE,
                show_column_dend = TRUE,
                cluster_column_slices = FALSE,
                row_title = "tumor-unique DMR (5329)", 
                column_title = "Samples",
                column_title_side = "bottom",
                column_names_gp = gpar(fontsize = 8),
                heatmap_legend_param = list(title = "Beta",
                                            title_position = "lefttop-rot",
                                            grid_height = unit(1, "cm"),
                                            grid_width = unit(0.5, "cm")))
  return(hm)
}

load("data/rdata/unique_lesions_filt.RData")
draw_hm(msetfem, sub_uniqueannot)

#### combine all markers into single signature and draw ROC ####

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
      colnames(meth_vals), median, na.rm=TRUE
    )
  
  #Get matrix
  score <- as.matrix(gr_dmr)[,-1]
  mean_betas <- colMedians(score, na.rm = TRUE)
}
mean_betas <- combine_marks(msetfem, sub_uniqueannot)

#set truth and draw ROC
truth <- ifelse(msetfem$tissue %in% c("adenoma","cancer"), 1,0)
pROC::roc(truth, mean_betas, direction = "<", plot = TRUE,
          main= "Combined selected markers, Luo 2014",
          percent=TRUE, print.auc=TRUE, print.thres = "best")


#calculate ROC for DMRs including age (les without seg, but with age)
load("data/rdata/les_with_age.RData")
mean_betas <- combine_marks(msetfem, sub_unique)

#set truth and draw ROC
truth <- ifelse(msetfem$tissue %in% c("adenoma","cancer"), 1,0)
pROC::roc(truth, mean_betas, direction = "<", plot = TRUE,
          main= "Combined selected markers with age, Luo 2014",
          percent=TRUE, print.auc=TRUE, print.thres = "best")
