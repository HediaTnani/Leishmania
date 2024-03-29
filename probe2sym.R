probe2sym <- function(sig, bmm, bmh) {
  require(biomaRt)
  genes <- getBM(attributes = c("affy_mogene_1_0_st_v1", "mgi_symbol","entrezgene_id","ensembl_gene_id"), filters = "affy_mogene_1_0_st_v1", values = rownames(sig), mart = bmm)
  m <- match(rownames(sig), genes$affy_mogene_1_0_st_v1)
  sig$affy <- genes[m, "affy"]
  sig$gene <- genes[m, "mgi_symbol"]
  sig$entrez_id<- genes[m, "entrezgene_id"]
  sig$ensembl_mm<- genes[m, "ensembl_gene_id"]
  x=sig$ensembl_mm
  genesV2 = getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = x , mart = bmm, attributesL = c("ensembl_gene_id"), martL = bmh,, uniqueRows=T)
  require(dplyr)
  require(readr)
  sig <- left_join(sig, genesV2, 
                   by = c("ensembl_mm" = "Gene.stable.ID"))
  sig$ensembl_hs=sig$Gene.stable.ID.1
  sig=dplyr::select(sig,-Gene.stable.ID.1)
  return(sig)
}
