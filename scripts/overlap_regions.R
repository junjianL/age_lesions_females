## Overlap Regions
# nov 26 2019

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ComplexHeatmap)
  library(ggplot2)
})


load("data/DMRs_age_cecum.RData")
load("data/DMRs_age_sig.RData")

#Get unique regions for each segment

findOverlaps(DMRsage_cecum, DMRsage_sig) #15,989
#cecum_unique <- subsetByOverlaps(DMRsage_cecum, DMRsage_sig, invert = TRUE)
#sig_unique <- subsetByOverlaps(DMRsage_sig, DMRsage_cecum, invert = TRUE)

#Filter and save files
DMRsage_cecum <- DMRsage_cecum[DMRsage_cecum$qval < 0.05] #674
DMRsage_sig <- DMRsage_sig[DMRsage_sig$qval < 0.05] #107

cecum_uniq <- subsetByOverlaps(DMRsage_cecum, DMRsage_sig, invert = TRUE) #632
sig_uniq <- subsetByOverlaps(DMRsage_sig, DMRsage_cecum, invert = TRUE)

df <- as.data.frame(cecum_uniq)
write.table(df, file = "data/DMRs_unique_cecum.txt", quote = FALSE, row.names = FALSE)

df <- as.data.frame(sig_uniq)
write.table(df, file = "data/DMRs_unique_sigmoid.txt", quote = FALSE, row.names = FALSE)

#Try ComplexHeatmap for upsets

DMRsage_cecum_hyper <- DMRsage_cecum[DMRsage_cecum$beta > 0] #656
DMRsage_cecum_hypo <- DMRsage_cecum[DMRsage_cecum$beta < 0] #18

DMRsage_sig_hyper <- DMRsage_sig[DMRsage_sig$beta > 0] #105
DMRsage_sig_hypo <- DMRsage_sig[DMRsage_sig$beta < 0] #2


makeUpsetTable <- function(x,y){
  hits <- findOverlaps(x, y)
  comm <- x[queryHits(hits)]
  
  uncomS <- y[-subjectHits(hits)]
  uncomA <- x[-queryHits(hits)]
  
  hyperm <- matrix(0, nrow = sum(length(comm), length(uncomA), length(uncomS)), ncol = 2)
  hyperm[0:length(comm),] <- 1
  hyperm[(length(comm)+1):(length(comm)+length(uncomA)),1] <- 1
  hyperm[(length(comm)+length(uncomA)+1):(length(comm)+length(uncomA)+length(uncomS)),2] <- 1
  
  
  colnames(hyperm)=c("Cecum (Old Vs Young)", "Sigmoid (Old Vs Young)")
  return(hyperm)
}

hyperM <- makeUpsetTable(DMRsage_cecum_hyper, DMRsage_sig_hyper)


#This works if i want to see number of bps overlapping
# peak_list_hyper <- list(Cecum = DMRsage_cecum_hyper, Sigmoid = DMRsage_sig_hyper)
# peak_list_hypo <- list(Cecum = DMRsage_cecum_hypo, Sigmoid = DMRsage_sig_hypo)

m = make_comb_mat(hyperM)
#m2 = make_comb_mat(peak_list_hypo) #no overlap

col_size = comb_size(m)
row_size = set_size(m)

pdf("figures/upset_hyper_cecumVssigmoid.pdf", width = 6, height = 5)
ht = UpSet(m, pt_size = unit(5, "mm"), lwd = 3,
           top_annotation = 
             HeatmapAnnotation("No. of DMRs\nhypermethylated" = anno_barplot(comb_size(m), 
                                                                             border = FALSE, 
                                                                             gp = gpar(fill = "#e1cf22"), 
                                                                             height = unit(7, "cm"),
                                                                             ylim = c(0, max(col_size)*1.1)
                               )),
           right_annotation = upset_right_annotation(m,
                                                     width = unit(4, "cm"),
                                                     ylim = c(0, max(row_size)*1.1)))
ht = draw(ht)

col_od = column_order(ht)
row_od = row_order(ht)

decorate_annotation("No. of DMRs\nhypermethylated", {
  grid.text(col_size[col_od], 
            seq_len(length(col_size)), 
            unit(col_size[col_od], "native") + unit(2, "mm"), 
            default.units = "native", just = "bottom",
            gp = gpar(fontsize = 8))
})
decorate_annotation("Set size", {
  grid.text(row_size[row_od], 
            unit(row_size[row_od], "native") + unit(2, "mm"), 
            rev(seq_len(length(row_size))), 
            default.units = "native", just = "bottom", rot = -90,
            gp = gpar(fontsize = 8))
})

dev.off()
