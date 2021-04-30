####################################################################
## Code to download data from dataset GSE113904, 
# and draw heatmaps
#
#
# August 4 2020
####################################################################
# Took some ideas from here https://f1000research.com/articles/5-1281/v3

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

# clean vars
mset$age <- as.integer(mset$`age:ch1`)
mset$age_group <- "mid"
mset$age_group[mset$age > 70] <- "old"
mset$age_group[mset$age <= 40] <- "young"
mset$age_group <- factor(mset$age_group, levels = c("young", "mid", "old"))
table(mset$age_group)

mset$location <- factor(mset$`anatomic location:ch1`)

mset$tissue <- factor(mset$`tissue type:ch1`)

msetfem <- mset[,mset$`gender:ch1` == "F"]
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


# save
saveRDS(msetfem, "data/public_data/luebeck.rds")

#mset <- readRDS("data/public_data/luebeck.rds")

#### Draw heatmaps ####

load("data/rdata/unique_lesions_filt.RData")

# get probes in selected markers
meth_vals <- as.matrix(getBeta(msetfem))

draw_hm <- function(obj, regions, splitby, ylab){
  gr <- rowRanges(obj)
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
  
  #print(gr_dmr[1:5,1:5])
  #Get matrix
  score <- as.matrix(gr_dmr)[,-1]
  print(dim(score))
  
  #return(score)
  
  #Colors
  col <- RColorBrewer::brewer.pal(n = 9, name = "YlGnBu")
  col_anot <- RColorBrewer::brewer.pal(n = 9, name = "Set1")
  purples <- RColorBrewer::brewer.pal(n = 3, name = "Purples")
  greens <- RColorBrewer::brewer.pal(n = 3, name = "Greens")
  pinks <- RColorBrewer::brewer.pal(n = 8, name = "RdPu")
  
  #hm colors
  col_fun <- circlize::colorRamp2(c(0,0.2,1), c(col[9], col[7], col_anot[6]))
  
  #annot colors
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
  
  #remove names
  rownames(score) <- colnames(score) <- NULL
  
  #Plot
  hm <- Heatmap(score, 
                use_raster = TRUE,
                na_col = "white",
                column_split = splitby,
                top_annotation = column_ha,
                col = col_fun,
                #row_km = 2, 
                clustering_distance_columns = "spearman",
                cluster_columns = TRUE,
                show_row_dend = FALSE,
                show_column_dend = TRUE,
                cluster_column_slices = FALSE,
                #left_annotation = row_ha, ## remove with 3 groups
                row_title = ylab,
                column_title = "Samples",
                column_title_side = "bottom",
                column_names_gp = gpar(fontsize = 8),
                heatmap_legend_param = list(title = "Beta",
                                            title_position = "lefttop-rot",
                                            grid_height = unit(1, "cm"),
                                            grid_width = unit(0.5, "cm")))
  return(hm)
}
draw_hm(msetfem, sub_uniqueannot, msetfem$tissue, "tumor-specific DMRs (5199)")

# draw age-specific DMRs
load("data/rdata/uniqueDMRs_age.RData")

draw_hm(msetfem, unique_age, msetfem$age_group, "age-specific DMRs (397)")

# draw age-associated DMRs
load("data/rdata/DMRs_age_final.RData")
age <- sort(DMRsage_annot[DMRsage_annot$qval <= 0.05]) #130 samps
draw_hm(msetfem, age, msetfem$age_group, "age-associated DMRs (2798)")

#remove mid-age to do hm
msetsub <- msetfem[,msetfem$age_group != "mid"]
meth_vals <- as.matrix(getBeta(msetsub))
draw_hm(msetsub, age, msetsub$age_group, "age-associated DMRs (2798)") #21


