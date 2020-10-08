###########################################################
## Upset plot from figure 3B 
# nov 26 2019
###########################################################

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ComplexHeatmap)
  library(ChIPpeakAnno)
})

cutoff <- 0.05

# get lesion regions
load("data/rdata/DMRs_lesions_3.RData")
lesion <- sort(DMRsles_annot[DMRsles_annot$qval <= cutoff & DMRsles_annot$beta > 0])

load("data/rdata/DMRs_lesions_SSA.RData")
ssa <- sort(DMRsles_annot[DMRsles_annot$qval <= cutoff & DMRsles_annot$beta > 0])

load("data/rdata/DMRs_lesions_cADN.RData")
cadn <- sort(DMRsles_annot[DMRsles_annot$qval <= cutoff & DMRsles_annot$beta > 0])

# get full age comparisons

load("data/rdata/DMRs_age_final.RData")
age_full <- DMRsage_annot[DMRsage_annot$qval <= cutoff]

#get full segment comparisons

load("data/rdata/DMRs_segm.RData")
seg_full <- DMRsage_annot[DMRsage_annot$qval <= cutoff]

# get overlap counts
# this returns input to venndiagram, I transform this to input
# for complexHeatmap

#function is based on the Interval tree algorithm
# over <- findOverlapsOfPeaks(lesion, ssa, cadn, age_full, seg_full, 
#                             connectedPeaks = "keepAll")
over <- findOverlapsOfPeaks(ssa, cadn, age_full, seg_full, 
                            connectedPeaks = "keepAll")

# mat <- data.frame(over$venn_cnt[,1:5])
# total_counts <- over$venn_cnt[,6]

mat <- data.frame(over$venn_cnt[,1:4])
total_counts <- over$venn_cnt[,5]


#Repeat rows `counts` number of times
fillin <- sapply(1:nrow(mat), function(u){
  if(total_counts[u] > 0) {
    reps <- unlist(rep(mat[u,], total_counts[u]))
    matrix(reps, ncol = ncol(mat), byrow = TRUE)
  } else NULL
})

fillin <- fillin[!sapply(fillin, is.null)]
matfull <- as.matrix(do.call(rbind,fillin))
colnames(matfull) <- colnames(mat)

#get obj for upset
m <- make_comb_mat(matfull)


#Plot
upset_withnumbers <- function(m, color, order_vec){
  col_size = comb_size(m)
  row_size = set_size(m)
  
  ht = UpSet(m, pt_size = unit(3, "mm"), lwd = 2, 
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
                                                       width = unit(3, "cm"),
                                                       ylim = c(0, max(row_size)*1.1),
                                                       gp = gpar(fill = "#636363")))
  ht = draw(ht)
  
  col_od = column_order(ht)
  row_od = row_order(ht)
  
  decorate_annotation("No. DMRs", {
    grid.text(col_size[col_od], 
              seq_len(length(col_size)), 
              unit(col_size[col_od], "native") + unit(1, "mm"), 
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

pdf("upsets_DMRs_3groups.pdf", width = 10, height = 6)

upset_withnumbers(m, "#756bb1", c("ssa", "cadn",
                                  "age_full","seg_full"))
dev.off()
