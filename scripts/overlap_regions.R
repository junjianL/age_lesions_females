## Overlap Regions
# nov 26 2019

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ComplexHeatmap)
  library(ggplot2)
})


#### Get regions that are uniquely in pre-lesions ####

cutoff <- 0.1

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

# get all combinations of age comparisons

load("data/rdata/DMRs_age_final.RData")
full <- sort(DMRsage_annot[DMRsage_annot$qval <= cutoff])
load("data/rdata/DMRs_age_sig_2.RData")
sig <- sort(DMRsage_sig_annot[DMRsage_sig_annot$qval <= cutoff])
load("data/rdata/DMRs_age_cecum_2.RData")
cec <- sort(DMRsage_cecum_annot[DMRsage_cecum_annot$qval <= cutoff])

age <- c(full, sig, cec)
age <- reduce(age)

#get all combinations of segment comparisons

load("data/rdata/DMRs_segm.RData")
full <- sort(DMRsage_annot[DMRsage_annot$qval <= cutoff])
load("data/rdata/DMRs_young.RData")
young <- sort(DMRsage_annot[DMRsage_annot$qval <= cutoff])
load("data/rdata/DMRs_old.RData")
old <- sort(DMRsage_annot[DMRsage_annot$qval <= cutoff])
rm(DMRsage_annot, DMRsage_cecum_annot, DMRsage_sig_annot)

seg <- c(full, young, old)
seg <- reduce(seg)

age_seg <- c(age, seg)
age_seg <- reduce(age_seg)

#get lesion regions without age or segment signal

unique_less <- subsetByOverlaps(les_source, age_seg, invert = TRUE)
unique_less
save(unique_less,file = sprintf("data/rdata/uniqueDMRs_lesion_%g.RData", cutoff))

df <- as.data.frame(unique_less)
write.table(df, file = sprintf("data/tables/uniqueDMRs_lesion_%g.txt",cutoff), 
            quote = FALSE, row.names = FALSE)


# merged
unique_less <- subsetByOverlaps(les, age_seg, invert = TRUE)
unique_less
save(unique_less,file = sprintf("data/rdata/uniqueDMRs_lesion_merged_%g.RData", cutoff))

df <- as.data.frame(unique_less)
write.table(df, file = sprintf("data/tables/uniqueDMRs_lesion_merged_%g.txt",cutoff), 
            quote = FALSE, row.names = FALSE)



#### Upset plots ####

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


