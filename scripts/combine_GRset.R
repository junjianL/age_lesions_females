## Combine all datasets
# August 4 2020

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(minfi)
  library(DMRcate)
})

# Import and combine
lueb <- readRDS("data/public_data/luebeck.rds")
wang <- readRDS("data/public_data/wang.rds")
diez <- readRDS("data/public_data/diez.rds")

over <- findOverlaps(lueb, wang)
lueb <- lueb[queryHits(over),]
wang <- wang[subjectHits(over),]
comb1 <- combine(lueb, wang)

over <- findOverlaps(comb1, diez)
diez <- diez[subjectHits(over),]
comb1 <- comb1[queryHits(over),]

assays(diez) <- assays(diez)[-2]
comb2 <- combine(comb1, diez)

# MDS plots

pal <- brewer.pal(8,"Dark2")

#by age
limma::plotMDS(getM(comb2), top=10000, gene.selection="common",
               col=pal[factor(comb2$age_group)], dim=c(1,3))
legend("top", legend=levels(factor(comb2$age_group)), text.col=pal,
       bg="white", cex=0.7)

#by dataset
limma::plotMDS(getM(comb2), top=10000, gene.selection="common",
               col=pal[factor(comb2$dataset)])
legend("top", legend=levels(factor(comb2$dataset)), text.col=pal,
       bg="white", cex=0.7)


# Run DMRcate
df <- colData(comb2)
designMatrix <- model.matrix(~ age_group + dataset, df)
myAnnotation <- cpg.annotate(datatype = "array", 
                             object = comb2, 
                             what = "M",
                             arraytype = "450K",
                             analysis.type="differential", 
                             design=designMatrix,
                             coef = 2)
str(myAnnotation)

DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
grDMRs <- extractRanges(DMRs, genome = "hg19")
save(grDMRs, file="data/rdata/450k_ageDMRs.RData")

# draw examples
groups <- pal[1:length(unique(comb2$age_group))]
names(groups) <- levels(factor(comb2$age_group))
cols <- groups[as.character(factor(comb2$age_group))]

# draw top DMR
DMR.plot(ranges=grDMRs, dmr=4, CpGs=comb2, phen.col=cols, what = "Beta",
         arraytype = "450K", pch=16, toscale=TRUE, plotmedians=TRUE, 
         genome="hg19")


# overlap DMRs with BS-seq age DMRs
load("data/rdata/DMRs_age_final.RData")
age <- DMRsage_annot[DMRsage_annot$qval <= 0.1]

sub <- subsetByOverlaps(age, grDMRs) #545, 636 with 0.1

# overlap DMRs with unique BS-seq lesion DMRs 
load("data/rdata/uniqueDMRs_lesion_merged_0.05.RData")
subbs <- subsetByOverlaps(unique_less, grDMRs)
subbs <- subsetByOverlaps(grDMRs,unique_less)
