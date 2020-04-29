probe2sym <- function(sig, bmm, bmh) {
  genes <- getBM(attributes = c("affy_mogene_1_0_st_v1", "mgi_symbol","entrezgene_id","ensembl_gene_id","hgnc_id"), filters = "affy_mogene_1_0_st_v1", values = rownames(sig), mart = bmm)
  m <- match(rownames(sig), genes$affy_mogene_1_0_st_v1)
  sig$gene <- genes[m, "mgi_symbol"]
  sig$entrez_id<- genes[m, "entrezgene_id"]
  x=sig$entrez_id
  genesV2 = getLDS(attributes = c("entrezgene_id"), filters = "entrezgene_id", values = x , mart = bmm, attributesL = c("hgnc_symbol"), martL = bmh)
  require(dplyr)
  require(readr)
  sig <- left_join(sig, genesV2, 
                   by = c("entrez_id" = "NCBI.gene.ID"))
  return(sig)
}
