## Public data
# August 4 2020

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(stringr)
  library(RColorBrewer)
  library(minfi)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
})

## DiezVillanueva data

#re-format metadata
metadata <- readr::read_delim("data/public_data/DiezVillanueva_2020/metadata.txt",
                              delim = "\t",
                              col_names = FALSE)
metadata <- t(metadata)[-1,]
metadata <- as.data.frame(metadata)
metadata <- metadata[,-6]
metadata$sample <- str_extract(metadata$V1, "[0-9A-Z]+_[A-Z]")

idxh <- grepl("healthy donor", metadata$V1)
hsamps <- metadata$sample[idxh]
metadata <- metadata[idxh,]

colnames(metadata) <- c("title", "tissue", "gender", "age", "location", "sample")

#clean var types
head(metadata)
metadata$age <- as.integer(gsub("age: ","",metadata$age))
metadata$age_group <- "none"
metadata$age_group[metadata$age > 70] <- "old"
metadata$age_group[metadata$age <= 40] <- "young"
table(metadata$age_group)
#only 1 sample under 40


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

#get healthy samps
pvals <- pvals[,hsamps]
betas <- betas[,hsamps]

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
                  diez$age_group != "none"]

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