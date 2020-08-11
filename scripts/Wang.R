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

metadata$age_group <- "none"
metadata$age_group[metadata$age > 70] <- "old"
metadata$age_group[metadata$age <= 40] <- "young"
table(metadata$age_group)

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

#subset females, age, low risk

wangfem <- wang[,wang$gender == "gender: F" & 
                  wang$age_group != "none" & 
                  wang$risk == "crc risk: Low"]

wangfem$dataset <- "Wang"
saveRDS(wangfem, "data/public_data/wang.rds")
#do some plots
pal <- brewer.pal(8,"Dark2")

#by age
limma::plotMDS(getM(wangfem), top=10000, gene.selection="common",
               col=pal[factor(wangfem$age_group)])
legend("top", legend=levels(factor(wangfem$age_group)), text.col=pal,
       bg="white", cex=0.7)

