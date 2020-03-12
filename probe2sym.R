probe2sym <- function(sig, bm) {
  genes <- getBM(attributes = c("affy_mogene_1_0_st_v1", "mgi_symbol"), filters = "affy_mogene_1_0_st_v1", values = rownames(sig), mart = bm)
  m <- match(rownames(sig), genes$affy_mogene_1_0_st_v1)
  sig$gene <- genes[m, "mgi_symbol"]
  return(sig)
}