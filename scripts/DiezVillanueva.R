####################################################################
## Code to download data from dataset GSE131013, 
# and draw figure 4B
#
# first see file `data/public_data/DiezVillanueva_2020/README.md`
# for file download
#
# August 4 2020
####################################################################

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
#diez <- readRDS("data/public_data/diez.rds")

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
source("helpers.R")

# get probes in selected markers
meth_vals <- as.matrix(getBeta(diezfem))

draw_hm <- function(obj, meth_vals, regions, splitby, ylab, metrics){
  
  score <- get_mat(obj, meth_vals, regions)[,-1]
  
  if(metric) mets <- get_metric(meth_vals, score) 
  
  #Colors
  col <- RColorBrewer::brewer.pal(n = 9, name = "YlGnBu")
  col_anot <- RColorBrewer::brewer.pal(n = 9, name = "Set1")
  purples <- RColorBrewer::brewer.pal(n = 3, name = "Purples")
  greens <- RColorBrewer::brewer.pal(n = 3, name = "Greens")
  pinks <- RColorBrewer::brewer.pal(n = 4, name = "RdPu")[c(1,3:4)]

  #hm colors
  col_fun <- circlize::colorRamp2(c(0,0.2,1), c(col[9], col[7], col_anot[6]))
  
  #row annot colors
  col_auc <- circlize::colorRamp2(c(0,50,100), c("#636363","#bdbdbd", "#dd1c77"))
  col_tpr <- circlize::colorRamp2(c(0,50,100), c("#636363","#bdbdbd", "#1c9099"))
  col_fpr <- circlize::colorRamp2(c(0,50,100), c("#636363","#bdbdbd", "#de2d26"))
  
  
  #col annot colors
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
  #row annot
  if(metrics){
  row_ha <- rowAnnotation("AUC" = mets$auc,
                          "TPR" = mets$sens,
                          "1-FPR" = mets$spec,
                          col = list("AUC" = col_auc,
                                     "TPR" = col_tpr,
                                     "1-FPR" = col_fpr))
  }
  
  #remove names
  rownames(score) <- colnames(score) <- NULL
  
  #Plot
  hm <- Heatmap(score, 
                use_raster = TRUE,
                na_col = "white",
                column_split = splitby,
                top_annotation = column_ha,
                col = col_fun,
                clustering_distance_columns = "spearman",
                cluster_columns = TRUE,
                cluster_rows = TRUE,
                show_row_dend = FALSE,
                show_column_dend = TRUE,
                cluster_column_slices = FALSE,
                #right_annotation = metrics, #turn on
                row_title = ylab, #5322
                column_title = "Samples",
                column_title_side = "bottom",
                heatmap_legend_param = list(title = "Beta",
                                            title_position = "lefttop-rot",
                                            grid_height = unit(1, "cm"),
                                            grid_width = unit(0.5, "cm")))
  return(hm)
}


#draw Figure 4.B
draw_hm(diezfem, sub_uniqueannot, diezfem$tissue,
        "Tumorigenesis-specific DMRs (5322)", TRUE) #329 samples

# draw heatmap for age-associated DMRs in normal tissue
dieznm <- diezfem[,diezfem$tissue != "Tumor"] #52 samples
dieznm$tissue <- droplevels(dieznm$tissue)

load("data/rdata/DMRs_age_final.RData")
age <- sort(DMRsage_annot[DMRsage_annot$qval <= 0.05]) #130 samps
meth_vals <- as.matrix(getBeta(dieznm))

draw_hm(dieznm, age, dieznm$age_group, "age_associated DMRs (2951)", FALSE)


# Add metrics to Supp.Table 1

metrics <- get_metric(meth_vals, beta_meds[,-1])

sub_uniqueannot$auc <- NA
sub_uniqueannot$auc[beta_meds[,1]] <- metrics[,1]


# get ROC curve for each DMR #

beta_meds <- get_mat(diezfem, meth_vals, sub_uniqueannot)

truth <- ifelse(grepl("_T",colnames(meth_vals)), 1,0)
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

point_df <- coords(test, "best")
point_df$phrase <- paste0(round(point_df$threshold,1)," (",
                           round(point_df$specificity,1),"%, ",
                          round(point_df$sensitivity,1),"%)")

# plot

ggplot() +
  scale_x_reverse() +
  geom_path(data = center_df, aes(fdr,tpr), size = 1.5) +
  geom_path(data = roc_df, aes(fdr,tpr, group = dmr), alpha = 0.005) +
  geom_point(data = point_df, aes(specificity,sensitivity), color = "blue", size = 5) +
  geom_label(data = point_df, aes(specificity-20,sensitivity-5, label = phrase)) +
  geom_label(aes(50,50, label = paste0("AUC: ",round(pROC::auc(test),1),"%"))) +
  theme_classic() +
  xlab("Specificity (%)") +
  ylab("Sensitivity (%)") +
  ggtitle("Combined selected markers, Diez-Villanueva 2020")

ggsave("Fig4B_ROC_diez.pdf", width = 5, height = 5)
