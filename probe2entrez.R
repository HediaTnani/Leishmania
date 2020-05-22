getGenes <- function(sig, bm) {
  genes <- getBM(attributes = c("affy_mogene_1_0_st_v1", "entrezgene_id"), 
                 filters = "affy_mogene_1_0_st_v1", 
                 values = rownames(sig), 
                 mart = bm)
  m <- match(rownames(sig), genes$affy_mogene_1_0_st_v1)
  sig$entrez_id <- genes[m, "entrezgene_id"]
  return(sig)
}
