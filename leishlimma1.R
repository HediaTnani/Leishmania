require(GEOquery)
require(oligo)
require(mogene10sttranscriptcluster.db)
require(RCurl)
require(foreign)
require(tidyverse)
require(affyPLM)
require(genefilter)
library(affycoretools)
setwd("~/")
celFiles <- list.celfiles("./gse31997data",full.names = TRUE)
url <- "https://raw.githubusercontent.com/HediaTnani/Leishmania/master/targets.csv"
targets <- getURL(url)                
my.targets <- read.AnnotatedDataFrame(textConnection(targets),header = TRUE, row.names = 1, sep=",") 
rawData <- read.celfiles(celFiles, phenoData = my.targets)
eset_rma <- oligo::rma(rawData)
eset <- getMainProbes(eset_rma)
filtered = nsFilter(eset, require.entrez=FALSE, remove.dupEntrez=FALSE)
eset_filtered <-filtered$eset
exp <- data.frame(exprs(eset_filtered))
expb <- exp[,c(1:18, 49:66)]
gse <- getGEO("GSE31997", AnnotGPL = TRUE)
pd <- pData(gse$GSE31997_series_matrix.txt.gz)
NIblb= as.character(pd$title[1:18])
NIC= as.character(pd$title[49:66])
pd.desNIb <- NIblb %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$")%>%
  as.data.frame() %>%
  unite("V2_V3", c("V2","V3"),sep = " " ) %>%
  dplyr::select(V2_V3,V4, V7, V8) %>%
  setNames(c("Genotype","Strain", "Time", "Replicate")) %>% 
  mutate(Genotype = gsub("Non Infected","NIBlb6",Genotype),
         Genotype = as.character(Genotype))

pd.desNIC <- NIC %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
  as.data.frame() %>% 
  unite("V2_V3", c("V2","V3"),sep = " " ) %>%
  dplyr::select(V2_V3,V4, V7, V8) %>%
  setNames(c("Genotype", "Strain", "Time", "Replicate")) %>% 
  mutate(Genotype = gsub("Non Infected","NIC57",Genotype),
         Genotype = as.character(Genotype))
pd.desb<- bind_rows(pd.desNIb,pd.desNIC)
rownames(pd.desb) <- pd$geo_accession[c(1:18, 49:66)]
Genotype = pd.desb$Genotype
eset_filtereds <- exprs(eset_filtered)[,c(1:18, 49:66)]
require(limma)
designMat<- model.matrix(~0+Genotype)
dim(designMat)
colnames(designMat) <- c("NIBlb6","NIC57")
cont.matrix <- makeContrasts(NIBlb6vsNIC57 = NIBlb6-NIC57,levels=designMat)
print(cont.matrix)
fit<-lmFit(expb, designMat)
fit.main<-contrasts.fit(fit, cont.matrix)
fit.main<-eBayes(fit.main)
class(fit.main)
topTab_NIBvsNIC <- topTable (fit.main, number=nrow(fit.main), coef="NIBlb6vsNIC57", adjust="fdr") 
topTab_NIBvsNIC_filt <- topTab_NIBvsNIC %>% filter(adj.P.Val < 0.05)
NIBvsNIC_filt_all<- topTab_NIBvsNIC_filt %>% filter(-1 > logFC  | logFC > 1)
nrow(NIBvsNIC_filt_all)
NIBvsNIC_filt_up <- topTab_NIBvsNIC_filt %>% filter( logFC > 1)
NIBvsNIC_filt_up
nrow(NIBvsNIC_filt_up)
NIBvsNIC_filt_down <- topTab_NIBvsNIC_filt %>% filter( logFC < -1)
nrow(NIBvsNIC_filt_down)
require(devtools)
source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2entrez.R")
require(biomaRt)
mart.mm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "mmusculus_gene_ensembl",
                   host="useast.ensembl.org")
NIBvsNIC_filt_all_annot   <- probe2entrez(NIBvsNIC_filt_all, mart.mm)
NIBvsNIC_filt_down_annot   <- probe2entrez(NIBvsNIC_filt_down, mart.mm)
NIBvsNIC_filt_up_annot   <- probe2entrez(NIBvsNIC_filt_up, mart.mm)

require(mogene10sttranscriptcluster.db)
geneSymbols <- select(mogene10sttranscriptcluster.db, rownames(topTab_NIBvsNIC_filt), c("SYMBOL"))
SYMBOLS<- geneSymbols$SYMBOL
volcanoplot(fit.main, coef=1, highlight=10, names=SYMBOLS, 
            main=paste("Differentially expressed genes", colnames(cont.matrix)[1], sep="\n"))
abline(v=c(-1,1))

library(org.Mm.eg.db)
library(clusterProfiler)
genelist=na.omit(as.character(NIBvsNIC_filt_all_annot$entrez_id))
eg = bitr(genelist, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Mm.eg.db")

ego2 <- enrichGO(gene         = eg$ENTREZID,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = "ENTREZID",
                 ont           = "ALL",
                 pAdjustMethod = "bonferroni",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
head(ego2,3)
dotplot(ego2, showCategory=3)

library(clusterProfiler)
search_kegg_organism('mmu', by='kegg_code')

sigGenes <- NIBvsNIC_filt_all_annot$entrez_id[ NIBvsNIC_filt_all_annot$P.Value < 0.05 & 
                                                 !is.na(NIBvsNIC_filt_all_annot$adj.P.Val) &
                                                 abs(NIBvsNIC_filt_all_annot$logFC) > 1 ]
sigGenes <- na.exclude(sigGenes)
kk <- enrichKEGG(gene = sigGenes, organism = 'mmu')
kk
head(kk, n=10)
library(pathview)
logFC <- NIBvsNIC_filt_all_annot$logFC
names(logFC) <- NIBvsNIC_filt_all_annot$entrez_id
pathview(gene.data = logFC, 
         pathway.id = "mmu04610", 
         species = "mmu", 
         limit = list(gene=5, cpd=1))

