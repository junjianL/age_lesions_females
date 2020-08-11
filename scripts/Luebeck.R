## Public data
# August 4 2020
# Took some ideas from here https://f1000research.com/articles/5-1281/v3

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(stringr)
  library(RColorBrewer)
  library(minfi)
})

## Luebeck data
#The 232 samples should all be healthy people

#### metadata ####
#re-format metadata
metadata <- readr::read_delim("data/public_data/Luebeck_2019/metadata.txt",
                              delim = "\t",
                              col_names = FALSE)
metadata <- t(metadata)[-1,]
metadata <- as.data.frame(metadata)

colnames(metadata) <- c("title", "tissue", "age", "gender", "location")

#clean var types
head(metadata)
metadata$sample <- str_extract(metadata$title, "SAMPLE_[0-9]+")
metadata$age <- as.integer(gsub("age: ","",metadata$age))

#enough females below 50?
table(metadata$age < 40 , metadata$gender)
#       gender: F gender: M
# FALSE       100        88
# TRUE         30        14


#### signal matrix ####

#import directly from GEO (this is basically an SE object)
mset <- getGenomicRatioSetFromGEO("GSE113904")

#subset females and clean vars
# also remove mid-age samples

mset$age <- as.integer(mset$`age:ch1`)
mset$age_group <- "none"
mset$age_group[mset$age > 70] <- "old"
mset$age_group[mset$age <= 40] <- "young"
table(mset$age_group)

mset$location <- mset$`anatomic location:ch1`

msetfem <- mset[,mset$`gender:ch1` == "F" & mset$age_group != "none"]
table(msetfem$age_group)

#do some plots
pal <- brewer.pal(8,"Dark2")

#by age
limma::plotMDS(getM(msetfem), top=10000, gene.selection="common",
        col=pal[factor(msetfem$age_group)])
legend("top", legend=levels(factor(msetfem$age_group)), text.col=pal,
       bg="white", cex=0.7)

#by location
limma::plotMDS(getM(msetfem), top=10000, gene.selection="common",
               col=pal[factor(msetfem$location)])
legend("top", legend=levels(factor(msetfem$location)), text.col=pal,
       bg="white", cex=0.7)

## remove rectum samples
msetfemsub <- msetfem[,!msetfem$location %in% c("rectum","right_colon")]

msetfemsub$dataset <- "Luebeck"
saveRDS(msetfemsub, "data/public_data/luebeck.rds")

#Run bumphunter wrapper in minfi
df <- colData(msetfem)
designMatrix <- model.matrix(~ age_group, df)

dmrs <- minfi::bumphunter(msetfem, design = designMatrix, 
                   cutoff = 0.1, B=1000, type="Beta")

dim(dmrs$table)

#load("data/public_data/Luebeck_2019/dmrs.RData")
