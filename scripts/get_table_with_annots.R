#!/usr/bin/env Rscript

#########################################################################################
# R function to add annotations from annotatr to a GRanges
#
# MBD-seq with 7 metastatis untreated, 5 metastasis treated, 6 CRC, 3 normal mucosa
# Collaboration with Mirco Menigatti and Giancarlo Marra
#
# Stephany Orjuela, February 2019
#########################################################################################

get_table_with_annots <- function(grreg){
  
  seqlevels(grreg) <- paste0("chr",seqlevels(grreg))
  annotscpg = c("hg19_cpg_islands", "hg19_cpg_shores",
                "hg19_cpg_shelves", "hg19_cpg_inter")
  
  annotations = build_annotations(genome = 'hg19', annotations = annotscpg)
  
  dm_annotated = annotate_regions(
    regions = grreg,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)
  
  df_dm_annotated <- data.frame(dm_annotated)
  
  over <- findOverlaps(grreg, dm_annotated)
  uniQ <- unique(queryHits(over))
  Q <- queryHits(over)
  S <- subjectHits(over)
  
  cpgAnnot <- sapply(uniQ, function(w){
    sub <- df_dm_annotated[S[Q == w],]
    paste(unique(sub$annot.id), collapse="/")
  })
  
  grreg$CpG <- "NA"
  grreg$CpG[uniQ] <- cpgAnnot
  
  #genes
  annotsgene <- c("hg19_genes_promoters", "hg19_genes_3UTRs", "hg19_genes_introns", 
                  "hg19_genes_exons", "hg19_genes_5UTRs", "hg19_genes_cds", "hg19_genes_intergenic",
                  "hg19_genes_1to5kb")
  
  annotations = build_annotations(genome = 'hg19', annotations = annotsgene)
  
  dm_annotated = annotate_regions(
    regions = grreg,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)
  
  df_dm_annotated <- data.frame(dm_annotated)
  
  over <- findOverlaps(grreg, dm_annotated)
  uniQ <- unique(queryHits(over))
  Q <- queryHits(over)
  S <- subjectHits(over)
  
  geneAnnot <- sapply(uniQ, function(w){
    sub <- df_dm_annotated[S[Q == w],]
    splited <- limma::strsplit2(sub$annot.id, ":")
    j <- paste(unique(paste(splited[,1], sub$annot.symbol, sep = ".")), collapse = "/")
    return(j)  
  })
  
  grreg$gene <- "NA"
  grreg$gene[uniQ] <- geneAnnot
  return(grreg)
}
