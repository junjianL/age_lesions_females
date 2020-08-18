## Upset plots 
# nov 26 2019

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ComplexHeatmap)
})

#Try ComplexHeatmap for upsets

cutoff <- 0.05

# get lesion regions
load("data/rdata/DMRs_lesions_3.RData")
lesion <- sort(DMRsles_annot[DMRsles_annot$qval <= cutoff & DMRsles_annot$beta > 0])

# get full age comparisons

load("data/rdata/DMRs_age_final.RData")
age_full <- DMRsage_annot[DMRsage_annot$qval <= cutoff & DMRsage_annot$beta > 0]

#get full segment comparisons

load("data/rdata/DMRs_segm.RData")
seg_full <- DMRsage_annot[DMRsage_annot$qval <= cutoff & DMRsage_annot$beta > 0]


#old function from old paper, too lazy to change
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

age_over <- makeUpsetTable(lesion, age_full, 
                         c("Lesion (Lesion Vs Normal)", "Age (Old Vs Young)"))

seg_over <- makeUpsetTable(lesion, seg_full, 
                        c("Lesion (Lesion Vs Normal)", "Location (Cecum Vs Sigmoid)"))


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

pdf("upsets_hyperDMRs.pdf")
m <- make_comb_mat(age_over)
upset_withnumbers(m, "#e1cf22", c("Lesion (Lesion Vs Normal)", "Age (Old Vs Young)"))

m <- make_comb_mat(seg_over)
upset_withnumbers(m, "#e1cf22", c("Lesion (Lesion Vs Normal)", "Location (Cecum Vs Sigmoid)")) 
dev.off()
#"#3b56d8" azul