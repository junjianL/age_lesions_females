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
})

## Wang data
#all samples are from left colon

#### metadata ####
#re-format metadata
metadata <- readr::read_delim("data/public_data/Wang_2020/metadata.txt",
                              delim = "\t",
                              col_names = FALSE)
metadata <- t(metadata)[-1,]
metadata <- as.data.frame(metadata)
colnames(metadata) <- c("sample", "ID", "source", "platform", 
                        "date", "age", "gender", "risk")

metadata$age <- as.integer(gsub("age: ","",metadata$age))

metadata$age_group <- "mid-age"
metadata$age_group[metadata$age > 70] <- "old"
metadata$age_group[metadata$age <= 40] <- "young"
metadata$age_group <- factor(metadata$age_group, levels = c("young","mid-age","old"))
table(metadata$age_group)

metadata$risk <- factor(gsub("crc risk: ","",metadata$risk), levels = c("Low","Medium","High"))

#### build GRSet ####
#Didnt work getting directly
#mset <- getGenomicRatioSetFromGEO("GSE132804", i = 1)
#Error in file(fname, "r") : invalid 'description' argument ??

# Import matrix parsed manually because i hate geo

betas <- read.table("data/public_data/Wang_2020/beta_mat.txt", 
                    header = TRUE, row.names = 1)

# get the 450k annotation data
ann450k <- getLocations(IlluminaHumanMethylation450kanno.ilmn12.hg19,
                        orderByLocation = TRUE)

ind2 <- match(rownames(betas),names(ann450k))

#build GRset
wang <- GenomicRatioSet(ann450k[ind2], 
                        Beta = betas,
                        annotation = c(array = "IlluminaHumanMethylation450k",
                                       annotation = "ilmn12.hg19"),
                        colData = metadata)

colnames(wang) <- colnames(betas)
wang$dataset <- "Wang"
saveRDS(wang, "data/public_data/wang.rds")

#subset females, age, low risk
wangfem <- wang[,wang$gender == "gender: F" & 
                  wang$age_group != "none" & 
                  wang$risk == "crc risk: Low"]

#do some plots
pal <- brewer.pal(8,"Dark2")

#by age
limma::plotMDS(getM(wangfem), top=10000, gene.selection="common",
               col=pal[factor(wangfem$age_group)])
legend("top", legend=levels(factor(wangfem$age_group)), text.col=pal,
       bg="white", cex=0.7)

#by sex
limma::plotMDS(getM(wangfem), top=10000, gene.selection="common",
               col=pal[factor(wangfem$gender)])
legend("top", legend=levels(factor(wangfem$gender)), text.col=pal,
       bg="white", cex=0.7)



#### Draw heatmap ####

wangfem <- wang[,wang$gender == "gender: F"]

load("data/rdata/unique_lesions_filt.RData")
#unique_all <- sub_unique[lengths(sub_unique$revmap) == 3] #2410
real <- overlapsAny(wangfem, sub_unique)
score <- data.frame(getBeta(wangfem)[real,])

  #Colors
col <- RColorBrewer::brewer.pal(n = 9, name = "YlGnBu")
col_anot <- RColorBrewer::brewer.pal(n = 9, name = "Set1")
purples <- RColorBrewer::brewer.pal(n = 3, name = "Purples")
oranges <- RColorBrewer::brewer.pal(n = 3, name = "Oranges")

  #hm colors
col_fun <- circlize::colorRamp2(c(0,0.2,1), c(col[9], col[7], col_anot[6]))
  
  #annot colors
col_age <- purples[1:nlevels(wangfem$age_group)]
names(col_age) <- levels(wangfem$age_group)
  
col_diag <- oranges[1:nlevels(wangfem$risk)]
names(col_diag) <- levels(wangfem$risk)
  
column_ha <- HeatmapAnnotation(Age = wangfem$age_group, 
                               Diagnosis = wangfem$risk,
                               col = list(Age = col_age,
                                            Diagnosis = col_diag), 
                               gp = gpar(col = "black"))

#remove names
rownames(score) <- colnames(score) <- NULL

  #Plot
Heatmap(score, 
        na_col = "white",
        column_split = wangfem$risk,
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

