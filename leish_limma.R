require(GEOquery)
require(oligo)
require(mogene10sttranscriptcluster.db)
require(RCurl)
require(foreign)
require(tidyverse)
require(affyPLM)
setwd("~/")
celFiles <- list.celfiles("./gse31997data",full.names = TRUE)
url <- "https://raw.githubusercontent.com/HediaTnani/Leishmania/master/targets.csv"
targets <- getURL(url)                
my.targets <- read.AnnotatedDataFrame(textConnection(targets),header = TRUE, row.names = 1, sep=",") 
rawData <- read.celfiles(celFiles, phenoData = my.targets)
eset_rma <- oligo::rma(rawData)
filtered = nsFilter(eset_rma, require.entrez=FALSE, remove.dupEntrez=FALSE)
eset_filtered <-filtered$eset
exp <- data.frame(exprs(eset_filtered))
expb <- exp[,4:33]
gse <- getGEO("GSE31997", AnnotGPL = TRUE)
pd <- pData(gse$GSE31997_series_matrix.txt.gz)
NIb= as.character(pd$title[4:18])
Ib= as.character(pd$title[19:33])
pd.desNIb <- NIb %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$")%>%
  as.data.frame() %>%
  unite("V2_V3", c("V2","V3"),sep = " " ) %>%
  dplyr::select(V2_V3,V7, V8) %>%
  setNames(c("Genotype","Time", "Replicate")) 
pd.desIb <- Ib %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>%
  as.data.frame() %>%
  dplyr::select(V8,V10, V11) %>%
  setNames(c("Genotype","Time", "Replicate")) %>% 
  mutate(Genotype = gsub("[[:punct:]]*", "", Genotype),
         Genotype = gsub("P", "Infected", Genotype),
         Genotype = as.character(Genotype))
pd.desb<- bind_rows(pd.desNIb,pd.desIb)
rownames(pd.desb) <- pd$geo_accession[4:33]
Genotype = pd.desb$Genotype
eset_filtereds <- exprs(eset_filtered)[,4:33]
require(limma)
designMat<- model.matrix(~0+Genotype)
dim(designMat)
colnames(designMat) <- c("NI","I")
cont.matrix <- makeContrasts(IvsNI = I-NI,levels=designMat)
print(cont.matrix)
fit<-lmFit(expb, designMat)
fit.main<-contrasts.fit(fit, cont.matrix)
fit.main<-eBayes(fit.main)
class(fit.main)
topTab_IvsNI <- topTable (fit.main, number=nrow(fit.main), coef="IvsNI", adjust="fdr") 
topTab_IvsNI_filt <- topTab_IvsNI %>% filter(adj.P.Val < 0.05)
IvsNI_filt_up <- topTab_IvsNI %>% filter( logFC > 1)
IvsNI_filt_up
nrow(IvsNI_filt_up)
IvsNI_filt_down <- topTab_IvsNI %>% filter( logFC < -1)
nrow(IvsNI_filt_down)
require(devtools)
source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2entrez.R")
require(biomaRt)
mart.mm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "mmusculus_gene_ensembl",
                   host = 'www.ensembl.org')
IvsNI_filt_down_annot   <- probe2entrez(IvsNI_filt_down, mart.mm)
IvsNI_filt_up_annot   <- probe2entrez(IvsNI_filt_up, mart.mm)

require(mogene10sttranscriptcluster.db)
geneSymbols <- select(mogene10sttranscriptcluster.db, rownames(topTab_IvsNI_filt), c("SYMBOL"))
SYMBOLS<- geneSymbols$SYMBOL
volcanoplot(fit.main, coef=1, highlight=10, names=SYMBOLS, 
            main=paste("Differentially expressed genes", colnames(cont.matrix)[1], sep="\n"))
abline(v=c(-1,1))
library(org.Mm.eg.db)
library(clusterProfiler)
genelist=na.omit(IvsNI_filt_down_annot$entrez_id)
ggo <- groupGO(gene     = as.character(genelist),
               OrgDb    = org.Mm.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)
head(ggo)
egoBP <- enrichGO(gene          = as.character(genelist),
                  OrgDb         = org.Mm.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
head(egoBP)
egoCC <- enrichGO(gene          = as.character(genelist),
                  OrgDb         = org.Mm.eg.db,
                  ont           = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
head(egoCC)
egoMF <- enrichGO(gene          = as.character(genelist),
                  OrgDb         = org.Mm.eg.db,
                  ont           = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
head(egoMF)
dotplot(egoBP, showCategory=30)
dotplot(egoCC, showCategory=30)
dotplot(egoMF, showCategory=10)


library("readxl")
PvsC_B = read_excel("~/Downloads/PvsC_B_v0.xlsx", sheet = "PvsC_B_v0")
PvsC_C = read_excel("~/Downloads/PvsC_C_v0.xlsx", sheet = "PvsC_C_v0")

m = PvsC_B$entrez_id %in% PvsC_C$entrez_id
PvsC_B[m,]
diff=setdiff(PvsC_B$entrez_id,PvsC_C$entrez_id)
diff
intersect=intersect(PvsC_B$entrez_id,PvsC_C$entrez_id)
diff1=na.omit(diff)
nrow(diff1)
intersect1=na.omit(intersect)
egoBP <- enrichGO(gene          = as.character(diff1),
                  OrgDb         = org.Mm.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
head(egoBP)
egoCC <- enrichGO(gene          = as.character(diff1),
                  OrgDb         = org.Mm.eg.db,
                  ont           = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
head(egoCC)
egoMF <- enrichGO(gene          = as.character(diff1),
                  OrgDb         = org.Mm.eg.db,
                  ont           = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
head(egoMF)
dotplot(egoBP, showCategory=30)
dotplot(egoCC, showCategory=30)
dotplot(egoMF, showCategory=30)
