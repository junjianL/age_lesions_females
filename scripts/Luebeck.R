## Public data
# August 4 2020

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(stringr)
})

## Luebeck data
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113904
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
table(metadata$age <= 50 , metadata$gender)
#       gender: F gender: M
# FALSE       100        88
# TRUE         30        14


#### signal matrix ####
library(minfi)

#import directly from GEO (this is basically an SE object)
mset <- getGenomicRatioSetFromGEO("GSE113904")

#subset females and clean vars

mset$age <- as.integer(mset$`age:ch1`)
mset$age_group <- as.factor(ifelse(mset$age > 50, "old", "young"))
table(mset$age_group)

mset$location <- mset$`anatomic location:ch1`

msetfem <- mset[,mset$`gender:ch1` == "F"]
table(msetfem$age_group)

#Run bumphunter wrapper in minfi
df <- colData(msetfem)
designMatrix <- model.matrix(~ age_group, df)

dmrs <- minfi::bumphunter(msetfem, design = designMatrix, 
                   cutoff = 0.1, B=1000, type="Beta")

dim(dmrs$table)


