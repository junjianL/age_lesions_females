####################################################################
## Code to download data from dataset TCGA - COAD
# and draw heatmaps
#
# July 30 2020
####################################################################


suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(stringr)
  library(RColorBrewer)
  library(curatedTCGAData)
  library(TCGAutils)
  library(ComplexHeatmap)
  library(dplyr)
  library(plyranges)
  
})

#### Use TCGA package to handle data

curatedTCGAData("COAD")

acc <- curatedTCGAData(
  diseaseCode = "COAD",
  assays = c(
    "Methylation_methyl450"
  ),
  dry.run = FALSE
)

sampleTables(acc)
#01  02  06  11   01: Primary Solid Tumor, 11: Solid Tissue Normal 
#293   1   1  38

getClinicalNames("COAD")

#choose females
#idx <- colData(acc)[["gender.x"]] == "female" | colData(acc)[["gender.y"]] == "FEMALE" 
#accfilt <- acc[,idx,]

#choose samples with ages
idx <- !is.na(colData(acc)[["years_to_birth"]])
accfilt <- acc[,idx,]

# Converting Assays to SummarizedExperiment
methcoad <- CpGtoRanges(accfilt)
rmethcoad <- experiments(methcoad)[[1]]

#add colData to RangedSCE
primaries <- unique(sampleMap(accfilt)$primary)
idx <- match(sampleMap(accfilt)$primary, primaries)

colData(rmethcoad) <- colData(accfilt)[idx,c("years_to_birth", "pathologic_stage","pathology_M_stage",
                                             "patient.gender")]
colData(rmethcoad)$ID <- sampleMap(accfilt)$colname
colData(rmethcoad)$tissue <- TCGAbiospec(sampleMap(accfilt)$colname)$sample_definition
rmethcoad <- rmethcoad[,!rmethcoad$tissue %in% c("Recurrent Solid Tumor","Metastatic")]

table(rmethcoad$tissue)

rmethcoad$tissue <- factor(rmethcoad$tissue,levels = c("Solid Tissue Normal", "Primary Solid Tumor"))

rmethcoad$age_group <- "mid-age"
rmethcoad$age_group[rmethcoad$years_to_birth > 70] <- "old"
rmethcoad$age_group[rmethcoad$years_to_birth <= 40] <- "young"
rmethcoad$age_group <- factor(rmethcoad$age_group, levels = c("young", "mid-age", "old"))
table(rmethcoad$age_group)

rmethcoad$patient.gender <- factor(rmethcoad$patient.gender)


#### Draw heatmap ####

source("scripts/helpers.R")

draw_hm <- function(obj, regions, meth_vals, ylab){
  
  score <- get_mat(obj, meth_vals, regions)[,-1]
  
  if(metric) mets <- get_metric(meth_vals, score)
  
  #Colors
  col <- RColorBrewer::brewer.pal(n = 9, name = "YlGnBu")
  col_anot <- RColorBrewer::brewer.pal(n = 9, name = "Set1")
  purples <- RColorBrewer::brewer.pal(n = 3, name = "Purples")
  greens <- RColorBrewer::brewer.pal(n = 4, name = "Greens")
  pinks <- RColorBrewer::brewer.pal(n = 3, name = "RdPu")
  
  #hm colors
  col_fun <- circlize::colorRamp2(c(0,0.2,1), c(col[9], col[7], col_anot[6]))
  
  #row annot colors
  col_auc <- circlize::colorRamp2(c(0,50,100), c("#636363","#bdbdbd", "#dd1c77"))
  col_tpr <- circlize::colorRamp2(c(0,50,100), c("#636363","#bdbdbd", "#1c9099"))
  col_fpr <- circlize::colorRamp2(c(0,50,100), c("#636363","#bdbdbd", "#de2d26"))
  
  #annot colors
  col_tis <- greens[1:nlevels(obj$tissue)]
  names(col_tis) <- levels(obj$tissue)
  
  col_age <- purples[1:nlevels(obj$age_group)]
  names(col_age) <- levels(obj$age_group)
  
  col_gen <- pinks[1:nlevels(obj$patient.gender)]
  names(col_gen) <- levels(obj$patient.gender)
  
  #column annot
  column_ha <- HeatmapAnnotation(
    Age = obj$age_group, 
    Tissue = obj$tissue,
    Gender = obj$patient.gender,
    col = list(
      Age = col_age,
      Tissue = col_tis,
      Gender = col_gen))
  
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
                show_column_dend = FALSE,
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
betas <- as.matrix(assay(rmethcoad, "counts"))
colnames(betas) <- rmethcoad$ID

#remove cpg with only Nas
idx <- rowSums(is.na(betas)) != ncol(betas)
betas <- betas[idx,]
rmethcoad <- rmethcoad[idx,]

draw_hm(rmethcoad, sub_uniqueannot, betas, "tumorigenesis-specific DMRs (5190)")

### draw hm for normal tissue, split by age
load("data/rdata/uniqueDMRs_age.RData")
rmethcoadnm <- rmethcoad[,rmethcoad$tissue != "Primary Solid Tumor"] 

#filter
betas <- as.matrix(assay(rmethcoadnm, "counts"))
colnames(betas) <- rmethcoadnm$ID

#remove cpg with only Nas
idx <- rowSums(is.na(betas)) != ncol(betas)
betas <- betas[idx,]
rmethcoadnm <- rmethcoadnm[idx,]

#plot
png("hmage_TCGA.png", width = 800, height = 700)
draw_hm(rmethcoadnm, unique_age, betas) #399 regions
dev.off()


#### combine all markers into single signature and draw ROC ####
seqlevels(rowRanges(rmethcoad)) <- paste0("chr",seqlevels(rowRanges(rmethcoad)))

beta_meds <- get_mat(rmethcoad, betas, sub_uniqueannot)

truth <- ifelse(rmethcoad$tissue == "Primary Solid Tumor", 1,0)
sens_mat <- apply(beta_meds, 1, function(u){
  rocc <- pROC::roc(truth, u[-1], direction = "<", plot = FALSE, percent = TRUE, quiet = TRUE)
  return(data.frame(fdr = rocc$specificities,
                    tpr = rocc$sensitivities, 
                    dmr = u[1]))
})


roc_df <- do.call(rbind, sens_mat)
roc_df$fdr <- jitter(roc_df$fdr, factor = 10, amount = 1)
roc_df$tpr <- jitter(roc_df$tpr, factor = 10, amount = 1)


# combine all markers into single signature and bolder ROC #
mean_betas <- colMedians(beta_meds[,-1], na.rm = TRUE)
test <- pROC::roc(truth, mean_betas, direction = "<", plot = TRUE,
                  percent=TRUE, print.auc=TRUE, print.thres = "best")

center_df <- data.frame(fdr = test$specificities,
                        tpr = test$sensitivities)

point_df <- pROC::coords(test, "best")
point_df$phrase <- paste0(round(point_df$threshold,1)," (",
                          round(point_df$specificity,1),"%, ",
                          round(point_df$sensitivity,1),"%)")

# plot

ggplot() +
  scale_x_reverse() +
  geom_path(data = center_df, aes(fdr,tpr), size = 1.5) +
  geom_path(data = roc_df, aes(fdr,tpr, group = dmr), alpha = 0.005) +
  geom_point(data = point_df, aes(specificity,sensitivity), color = "blue", size = 5) +
  geom_label(data = point_df, aes(specificity-20,sensitivity-4, label = phrase)) +
  geom_label(aes(50,50, label = paste0("AUC: ",round(pROC::auc(test),1),"%"))) +
  theme_classic() +
  xlab("Specificity (%)") +
  ylab("Sensitivity (%)") +
  ggtitle("Combined selected markers, TCGA")

ggsave("ROC_TCGA.pdf", width = 5, height = 5)




