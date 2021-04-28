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

#choose adenocarcinomas
# idx <- colData(accfilt)[["histological_type.x"]] == "colon adenocarcinoma" | 
#   colData(accfilt)[["histological_type.y"]] == "Colon Adenocarcinoma" 
# accfilt <- accfilt[,idx,]

#split and choose samples with paired normal tissue
#(first way i figured how to do this)

#grab barcodes with sample type 11
# idx <- grepl("^11", TCGAbarcode(sampleMap(accfilt)$colname,
#                          participant = FALSE, sample = TRUE))

#grab sampleID from the above
# normsamps <- sampleMap(accfilt)$primary[idx]

#find barcodes containing these sampleIDs
# idx <- grepl(paste(normsamps,collapse="|"), sampleMap(accfilt)$colname)

#these are the patients with paired samples
# final <- unique(sampleMap(accfilt)[idx,"primary"])
# accfilt <- accfilt[,final,]

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

draw_hm <- function(obj, regions, meth_vals){
  gr <- rowRanges(obj)
  seqlevels(gr) <- paste0("chr",seqlevels(gr))
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
  
  gr_dmr
  #Get matrix
  madr <- rowVars(as.matrix(gr_dmr)[,-1])
  o <- order(madr, decreasing = TRUE)
  score <- as.matrix(gr_dmr)[o,-1]
  print(dim(score))
  
  #return(score)
  
  #Colors
  col <- RColorBrewer::brewer.pal(n = 9, name = "YlGnBu")
  col_anot <- RColorBrewer::brewer.pal(n = 9, name = "Set1")
  purples <- RColorBrewer::brewer.pal(n = 3, name = "Purples")
  greens <- RColorBrewer::brewer.pal(n = 4, name = "Greens")
  pinks <- RColorBrewer::brewer.pal(n = 3, name = "RdPu")
  
  #hm colors
  col_fun <- circlize::colorRamp2(c(0,0.2,1), c(col[9], col[7], col_anot[6]))
  
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
  
  #remove names
  rownames(score) <- colnames(score) <- NULL
  
  #Plot
  hm <- Heatmap(score, 
                use_raster = TRUE,
                na_col = "white",
                column_split = obj$age_group,
                top_annotation = column_ha,
                col = col_fun,
                clustering_distance_columns = "spearman",
                cluster_columns = TRUE,
                show_row_dend = FALSE,
                show_column_dend = FALSE,
                cluster_column_slices = FALSE,
                row_title = "age-unique DMR (399)", 
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

png("hm_TCGA.png", width = 800, height = 700)
draw_hm(rmethcoad, sub_uniqueannot, betas)
dev.off()

#draw hm for normal tissue, split by age
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

combine_marks <- function(obj, regions, meth_vals){
  gr <- rowRanges(obj)
  seqlevels(gr) <- paste0("chr",seqlevels(gr))
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
mean_betas <- combine_marks(rmethcoad, sub_uniqueannot, betas)

#set truth and draw ROC
truth <- ifelse(rmethcoad$tissue == "Primary Solid Tumor", 1,0)
pROC::roc(truth, mean_betas, direction = "<", plot = TRUE,
          main= "Combined selected markers, TCGA",
          percent=TRUE, print.auc=TRUE, print.thres = "best")


