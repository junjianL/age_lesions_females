## Overlap Regions
# nov 26 2019

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ComplexHeatmap)
  library(ggplot2)
})


#### Get unique regions for each segment ####
load("data/DMRs_age_cecum_2.RData")
load("data/DMRs_age_sig_2.RData")

#Filter

DMRsage_sig <- DMRsage_sig_annot[DMRsage_sig_annot$qval < 0.05] #112
DMRsage_cecum <- DMRsage_cecum_annot[DMRsage_cecum_annot$qval < 0.05] #697

#Filter and save files
cecum_uniq <- subsetByOverlaps(DMRsage_cecum, DMRsage_sig, invert = TRUE) #632, 651
sig_uniq <- subsetByOverlaps(DMRsage_sig, DMRsage_cecum, invert = TRUE) #66

df <- as.data.frame(cecum_uniq)
write.table(df, file = "data/DMRs_unique_cecum.txt", quote = FALSE, row.names = FALSE)

df <- as.data.frame(sig_uniq)
write.table(df, file = "data/DMRs_unique_sigmoid.txt", quote = FALSE, row.names = FALSE)

#Try ComplexHeatmap for upsets

DMRsage_cecum_hyper <- DMRsage_cecum[DMRsage_cecum$beta > 0] #674
DMRsage_cecum_hypo <- DMRsage_cecum[DMRsage_cecum$beta < 0] #23 

DMRsage_sig_hyper <- DMRsage_sig[DMRsage_sig$beta > 0] #110
DMRsage_sig_hypo <- DMRsage_sig[DMRsage_sig$beta < 0] #2


makeUpsetTable <- function(x,y, comp){
  hits <- findOverlaps(x, y)
  comm <- x[queryHits(hits)]
  
  uncomS <- y[-subjectHits(hits)]
  uncomA <- x[-queryHits(hits)]
  
  hyperm <- matrix(0, nrow = sum(length(comm), length(uncomA), length(uncomS)), ncol = 2)
  hyperm[0:length(comm),] <- 1
  hyperm[(length(comm)+1):(length(comm)+length(uncomA)),1] <- 1
  hyperm[(length(comm)+length(uncomA)+1):(length(comm)+length(uncomA)+length(uncomS)),2] <- 1
  
  
  #colnames(hyperm)=c("Segments (Old Vs Young)","Lesions (Vs Normal)")
  colnames(hyperm)=comp
  return(hyperm)
}

hyperM <- makeUpsetTable(DMRsage_cecum_hyper, DMRsage_sig_hyper, 
                         c("Cecum (Old Vs Young)", "Sigmoid (Old Vs Young)"))

hypoM <- makeUpsetTable(DMRsage_cecum_hypo, DMRsage_sig_hypo, 
                         c("Cecum (Old Vs Young)", "Sigmoid (Old Vs Young)"))


#This works if i want to see number of bps overlapping
# peak_list_hyper <- list(Cecum = DMRsage_cecum_hyper, Sigmoid = DMRsage_sig_hyper)
# peak_list_hypo <- list(Cecum = DMRsage_cecum_hypo, Sigmoid = DMRsage_sig_hypo)

#m2 = make_comb_mat(peak_list_hypo) #no overlap

#Plot
upset_withnumbers <- function(m, color, order_vec){
  col_size = comb_size(m)
  row_size = set_size(m)
  
  ht = UpSet(m, pt_size = unit(5, "mm"), lwd = 3, 
             set_order = order_vec,
             top_annotation = 
               HeatmapAnnotation("No. DMRs" = 
                                   anno_barplot(comb_size(m), 
                                                         border = FALSE, 
                                                         gp = gpar(fill = color), 
                                                         height = unit(7, "cm"),
                                                         ylim = c(0, max(col_size)*1.1)
               )),
             right_annotation = upset_right_annotation(m,
                                                       width = unit(4, "cm"),
                                                       ylim = c(0, max(row_size)*1.1)))
  ht = draw(ht)
  
  col_od = column_order(ht)
  row_od = row_order(ht)
  
  decorate_annotation("No. DMRs", {
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
}
m <- make_comb_mat(hyperM)
pdf("figures/upset_cecumVssigmoid_pval.pdf", width = 6, height = 5)
upset_withnumbers(m, "#e1cf22", c("Cecum (Old Vs Young)", "Sigmoid (Old Vs Young)"))
m <- make_comb_mat(hypoM)
upset_withnumbers(m, "#3b56d8", c("Cecum (Old Vs Young)", "Sigmoid (Old Vs Young)")) 
dev.off()


#### Overlap with lesion DMRs ####

load("data/DMRs_lesions_2.RData")
load("data/DMRs_age_2.RData")
DMRsage <- DMRsage_annot[DMRsage_annot$qval < 0.05]
DMRsles <- DMRsles_annot[DMRsles_annot$qval < 0.07]
length(DMRsles)
#upset test, all 4 granges, sets are base pairs
peak_list <- list(Cecum = DMRsage_cecum, 
                  Sigmoid = DMRsage_sig, 
                  Segments = DMRsage,
                  Lesions = DMRsles)

m2 = make_comb_mat(peak_list)

pdf("figures/upset_full_bps.pdf", width = 8, height = 6)
upset_withnumbers(m2, "black", c("Lesions","Segments","Cecum","Sigmoid"))
dev.off()

#lesions and age regions
mat <- makeUpsetTable(DMRsles,DMRsage, c("Lesions (Vs Normal)", "Segments (Old Vs Young)"))
m <- make_comb_mat(mat)
pdf("figures/upset_lesionsVssegments.pdf",width = 6, height = 5)
upset_withnumbers(m, "black",c("Lesions (Vs Normal)", "Segments (Old Vs Young)"))
dev.off()


#split by hyper and hypo
DMRsage_hyper <- DMRsage[DMRsage$beta > 0] 
DMRsage_hypo <- DMRsage[DMRsage$beta < 0] 

DMRsles_hyper <- DMRsles[DMRsles$beta > 0]
DMRsles_hypo <- DMRsles[DMRsles$beta < 0] 

mat_hyper <- makeUpsetTable(DMRsles_hyper, DMRsage_hyper,
                            c("Lesions (Vs Normal)", "Segments (Old Vs Young)"))
mat_hypo <- makeUpsetTable(DMRsles_hypo, DMRsage_hypo,
                           c("Lesions (Vs Normal)", "Segments (Old Vs Young)"))

m <- make_comb_mat(mat_hyper)
m0 <- make_comb_mat(mat_hypo)
pdf("figures/upset_lesionsVssegments_hyper_hypo.pdf",width = 6, height = 5)
upset_withnumbers(m, "#e1cf22", c("Lesions (Vs Normal)", "Segments (Old Vs Young)"))
upset_withnumbers(m0, "#3b56d8", c("Lesions (Vs Normal)", "Segments (Old Vs Young)"))
dev.off()

#Filter and save files
age_uniq <- subsetByOverlaps(DMRsage, DMRsles, invert = TRUE) 
lesions_uniq <- subsetByOverlaps(DMRsles, DMRsage, invert = TRUE) 

df <- as.data.frame(age_uniq)
write.table(df, file = "data/DMRs_age_not_in_lesions.txt", quote = FALSE, row.names = FALSE)

df <- as.data.frame(lesions_uniq)
write.table(df, file = "data/DMRs_lesions_not_in_age.txt", quote = FALSE, row.names = FALSE)
