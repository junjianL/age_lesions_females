####################################################################
## Helper functions for Figure 4
#
# May 7 2021
####################################################################


### summarize methylation across all sites that overlap a region/probe/DMR
get_mat <- function(obj, meth_vals, regions){
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
  
  #Get matrix
  score <- as.matrix(gr_dmr)
  print(dim(score))
  return(score)
  
}

### get ROC metrics for each row in a matrix of meth vals
get_metric <- function(meth_vals, score, truth){
  #Get AUC, TPR, 1-FPR
  
  auc <- apply(score, 1, function(u){
    rocc <- pROC::roc(truth, u, direction = "<", plot = FALSE, percent = TRUE, quiet = TRUE)
    pROC::auc(rocc)
  })
  
  sens <- apply(score, 1, function(u){
    rocc <- pROC::roc(truth, u, direction = "<", plot = FALSE, percent = TRUE, quiet = TRUE)
    thresh <- pROC::coords(rocc, "best", transpose = TRUE)
    thresh[[3]]
  })
  
  spec <- apply(score, 1, function(u){
    rocc <- pROC::roc(truth, u, direction = "<", plot = FALSE, percent = TRUE, quiet = TRUE)
    thresh <- pROC::coords(rocc, "best", transpose = TRUE)
    thresh[[2]]
  })
  
  return(cbind(auc,sens,spec))
}

