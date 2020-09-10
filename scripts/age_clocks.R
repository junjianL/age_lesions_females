## Age and clocks
# Sep 9 2020

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  
})


cutoff <- 0.05

# get all combinations of lesion comparisons

load("data/rdata/DMRs_lesions_3.RData")
full <- sort(DMRsles_annot[DMRsles_annot$qval <= cutoff & DMRsles_annot$beta > 0])
full$comparison <- "SSA+cADN"
load("data/rdata/DMRs_lesions_SSA.RData")
ssa <- sort(DMRsles_annot[DMRsles_annot$qval <= cutoff & DMRsles_annot$beta > 0])
ssa$comparison <- "SSA"
load("data/rdata/DMRs_lesions_cADN.RData")
cadn <- sort(DMRsles_annot[DMRsles_annot$qval <= cutoff & DMRsles_annot$beta > 0])
cadn$comparison <- "cADN"
rm(DMRsles_annot)

les_source <- c(full,ssa,cadn)
o <- order(les_source$qval)
les_source <- les_source[o]
les <- reduce(les_source, with.revmap=TRUE)

# get age comparison

load("data/rdata/DMRs_age_final.RData")
age <- sort(DMRsage_annot[DMRsage_annot$qval <= cutoff])


#get lesion regions without age or segment signal
# merged
unique_less <- subsetByOverlaps(les, age_seg, invert = TRUE)
unique_less
save(unique_less,file = sprintf("data/rdata/uniqueDMRs_lesion_merged_%g.RData", cutoff))

#### compare to epiTOC CpGs ####
load("data/public_data/13059_2016_1064_MOESM5_ESM.rdata")
ann450k <- getLocations(IlluminaHumanMethylation450kanno.ilmn12.hg19,
                        orderByLocation = TRUE)

ind2 <- match(epiTOCcpgs.v,names(ann450k))
epitoc <- ann450k[ind2]

subsetByOverlaps(epitoc, les) #255
subsetByOverlaps(epitoc, age) #63

#### compate to hovarth CpGs ####

probes <- read.csv("data/public_data/13059_2013_3156_MOESM3_ESM.csv", header = FALSE)
hovarth <- probes$V1[-(1:4)]

ind2 <- match(hovarth,names(ann450k))
hov <- ann450k[ind2]

subsetByOverlaps(hov, les) #56
subsetByOverlaps(hov, age) #9
