## combine all samples
# dec 16 2019

suppressPackageStartupMessages({
  library(dmrseq)
  library(ggplot2)
})

bsList <- list(bismarkBSseq, bismarkBSseq_les)
bsCombined <- combineList(bsList)

#Filter
loci.idx <- DelayedMatrixStats::rowSums2(getCoverage(bsCombined, type="Cov") >= 10 ) >= 38
bsCombined <- bsCombined[loci.idx,]

save(bsCombined, file = "data/bsseqCombined.RData")

#Get Cov and meth
cov <- getCoverage(bsCombined, type = "Cov")
meth <- getCoverage(bsCombined, type = "M")
meth_vals <- meth /cov
colnames(meth_vals) <- colData(bsCombined)$names
colData(bsCombined)$state <- ifelse(is.na(colData(bsCombined)$state), 
                                    "Normal", colData(bsCombined)$state)

idx <- colData(bsCombined)$state == "Normal"

#Transform
methsTR <- asin(2*meth_vals-1)[,idx]

#Make MDS
mds_meth <- limma::plotMDS(methsTR, top = 10000, plot = FALSE)$cmdscale.out

colData(bsCombined)$lesion <- ifelse(colData(bsCombined)$segment %in% c("cecum", "sigmoid"), 
                                     paste0(colData(bsCombined)$segment, "_", colData(bsCombined)$age_group), 
                                     colData(bsCombined)$lesion)

#Plot
df <- data.frame(dim1 = mds_meth[,1], dim2 = mds_meth[,2],
                 names = colData(bsCombined)$patient[idx], treat = colData(bsCombined)$lesion[idx])

myColor <- RColorBrewer::brewer.pal(9, "Set1")
ggplot()+
  geom_point(data = df, mapping = aes_(x=~dim1, y=~dim2, color=~treat), size = 7) +
  geom_text(data = df, mapping = aes_(x=~dim1, y=~dim2-0.01, label=~names), size = 3) +
  scale_color_manual(values = myColor) +
  theme_bw()
ggsave("figures/MDS_allsamps.png")