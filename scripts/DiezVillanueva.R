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

# choose females and age
diezfem <- diez[,diez$gender == "gender: Female" & 
                  diez$age_group != "mid-age"]

# remove low detection pvalues across all samples
idx <- rowSums(assay(diezfem, "pvalue") < 0.05) == 6
diezfem <- diezfem[idx,]

diezfem$dataset <- "DiezVillanueva"
saveRDS(diezfem, "data/public_data/diez.rds")

#do some plots
pal <- brewer.pal(8,"Dark2")

#by age
limma::plotMDS(getM(diezfem), top=10000, gene.selection="common",
               col=pal[factor(diezfem$age_group)])
legend("top", legend=levels(factor(diezfem$age_group)), text.col=pal,
       bg="white", cex=0.7)

#by location
limma::plotMDS(getM(diezfem), top=10000, gene.selection="common",
               col=pal[factor(diezfem$location)])
legend("top", legend=levels(factor(diezfem$location)), text.col=pal,
       bg="white", cex=0.7)

#### draw heatmap with selected markers ####

load("data/rdata/unique_lesions_filt.RData")
#unique_all <- sub_unique[lengths(sub_unique$revmap) == 3] #2410

diezfem <- diez[,diez$gender == "gender: Female"]
dim(diezfem)

# remove low detection pvalues across all samples
idx <- rowSums(assay(diezfem, "pvalue") < 0.05) == 78
diezfem <- diezfem[idx,]
dim(diezfem) #481,955

# get probes in selected markers
real <- overlapsAny(diezfem, sub_unique)
score <- data.frame(getBeta(diezfem)[real,])

#Colors
col <- RColorBrewer::brewer.pal(n = 9, name = "YlGnBu")
col_anot <- RColorBrewer::brewer.pal(n = 9, name = "Set1")
purples <- RColorBrewer::brewer.pal(n = 3, name = "Purples")
greens <- RColorBrewer::brewer.pal(n = 3, name = "Greens")
pinks <- RColorBrewer::brewer.pal(n = 4, name = "RdPu")[c(1,3:4)]

#hm colors
col_fun <- circlize::colorRamp2(c(0,0.2,1), c(col[9], col[7], col_anot[6]))

#annot colors
col_age <- purples[1:nlevels(diezfem$age_group)]
names(col_age) <- levels(diezfem$age_group)

col_tis <- greens[1:nlevels(diezfem$tissue)]
names(col_tis) <- levels(diezfem$tissue)

col_seg <- pinks[1:nlevels(diezfem$location)]
names(col_seg) <- levels(diezfem$location)

#column annot
column_ha <- HeatmapAnnotation(Age = diezfem$age_group, 
                               Tissue = diezfem$tissue,
                               Location = diezfem$location,
                               col = list(Age = col_age,
                                          Tissue = col_tis,
                                          Location = col_seg), 
                               gp = gpar(col = "black"))

#remove names
rownames(score) <- colnames(score) <- NULL

#Plot
Heatmap(score, 
        na_col = "white",
        column_split = diezfem$tissue,
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

