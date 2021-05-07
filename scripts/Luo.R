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

#### Draw heatmap Figure 4.A ####
meth_vals <- as.matrix(getBeta(msetfem))

source("helpers.R")

draw_hm <- function(obj, regions, meth_vals, ylab){
  score <- get_mat(obj, meth_vals, regions)[,-1]
  
  truth <- ifelse(obj$tissue %in% c("adenoma","cancer"), 1,0)
  
  if(metric) mets <- get_metric(meth_vals, score, truth) 
  
  #Colors
  col <- RColorBrewer::brewer.pal(n = 9, name = "YlGnBu")
  col_anot <- RColorBrewer::brewer.pal(n = 9, name = "Set1")
  greens <- RColorBrewer::brewer.pal(n = 4, name = "Greens")
  pinks <- RColorBrewer::brewer.pal(n = 7, name = "RdPu")
  
  #hm colors
  col_fun <- circlize::colorRamp2(c(0,0.2,1), c(col[9], col[7], col_anot[6]))
  
  #row annot colors
  col_auc <- circlize::colorRamp2(c(0,50,100), c("#636363","#bdbdbd", "#dd1c77"))
  col_tpr <- circlize::colorRamp2(c(0,50,100), c("#636363","#bdbdbd", "#1c9099"))
  col_fpr <- circlize::colorRamp2(c(0,50,100), c("#636363","#bdbdbd", "#de2d26"))
  
  #annot colors
  col_tis <- greens[1:nlevels(obj$tissue)]
  names(col_tis) <- levels(obj$tissue)
  
  col_seg <- pinks[1:nlevels(obj$location)]
  names(col_seg) <- levels(obj$location)
  
  #column annot
  column_ha <- HeatmapAnnotation( 
                                 Tissue = obj$tissue,
                                 Location = obj$location,
                                 col = list(
                                            Tissue = col_tis,
                                            Location = col_seg), 
                                 gp = gpar(col = "black"))
  #row annot
  row_ha <- rowAnnotation("AUC" = mets$auc,
                          "TPR" = mets$sens,
                          "1-FPR" = mets$spec,
                          col = list("AUC" = col_auc,
                                     "TPR" = col_tpr,
                                     "1-FPR" = col_fpr))
  
  #remove names
  rownames(score) <- colnames(score) <- NULL
  
  #Plot
  hm <- Heatmap(score, 
                use_raster = TRUE,
                na_col = "white",
                column_split = obj$tissue,
                top_annotation = column_ha,
                right_annotation = row_ha,
                col = col_fun,
                clustering_distance_columns = "spearman",
                cluster_columns = TRUE,
                show_row_dend = FALSE,
                show_column_dend = TRUE,
                cluster_column_slices = FALSE,
                row_title = ylab, 
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
draw_hm(msetfem, sub_uniqueannot, meth_vals, "Tumorigenesis-specific DMRs (5329)")

#### Get ROC curve for each DMR ####


multi_roc <- function(idx, msetfem, meth_vals, sub_uniqueannot, tissue){
  
  beta_meds <- get_mat(msetfem[,idx], meth_vals[,idx], sub_uniqueannot)
  truth <- ifelse(msetfem$tissue[idx] == tissue, 1,0)
  
  sens_mat <- apply(beta_meds, 1, function(u){
    rocc <- pROC::roc(truth, u[-1], direction = "<", plot = FALSE, percent = TRUE, quiet = TRUE)
    return(data.frame(fdr = rocc$specificities,
                      tpr = rocc$sensitivities, 
                      dmr = u[1]))
  })
  
  
  roc_df <- do.call(rbind, sens_mat)
  roc_df$fdr <- jitter(roc_df$fdr, factor = 10, amount = 1)
  roc_df$tpr <- jitter(roc_df$tpr, factor = 10, amount = 1)
  
  ### combine all markers into single signature and draw ROC ###
  
  mean_betas <- colMedians(beta_meds[,-1], na.rm = TRUE)
  
  test <- pROC::roc(truth, mean_betas, direction = "<", plot = FALSE,
                    percent=TRUE, print.auc=TRUE, print.thres = "best")
  
  center_df <- data.frame(fdr = test$specificities,
                          tpr = test$sensitivities)
  
  point_df <- pROC::coords(test, "best")
  point_df$phrase <- paste0(round(point_df$threshold,1)," (",
                            round(point_df$specificity,1),"%, ",
                            round(point_df$sensitivity,1),"%)")
  
  # plot
  
  p <- ggplot() +
    scale_x_reverse() +
    geom_path(data = center_df, aes(fdr,tpr), size = 1.5) +
    geom_path(data = roc_df, aes(fdr,tpr, group = dmr), alpha = 0.005) +
    geom_point(data = point_df, aes(specificity,sensitivity), color = "blue", size = 5) +
    geom_label(data = point_df, aes(specificity-20,sensitivity-4, label = phrase)) +
    geom_label(aes(50,50, label = paste0("AUC: ",round(pROC::auc(test),1),"%"))) +
    theme_classic() +
    xlab("Specificity (%)") +
    ylab("Sensitivity (%)") +
    ggtitle(sprintf("Combined selected markers, Luo 2014 (%s)",tissue))
  
  return(p)
}

#adenomas
pdf("Fig4A_ROCs_luo.pdf", width = 5, height = 5)
aden <- msetfem$tissue != "cancer"
multi_roc(aden, msetfem, meth_vals, sub_uniqueannot, "adenoma")

#cancer
crc <- msetfem$tissue != "adenoma"
multi_roc(crc, msetfem, meth_vals, sub_uniqueannot, "cancer")
dev.off()


