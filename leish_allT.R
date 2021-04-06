# All pairwise comparison in limma 
#### T0 ####
require(GEOquery)
require(oligo)
require(mogene10sttranscriptcluster.db)
require(RCurl)
require(foreign)
require(tidyverse)
require(affyPLM)
require(genefilter)
require(AnnotationDbi)
library(affycoretools)
setwd("~/")
#install.packages("http://mbni.org/customcdf/25.0.0/entrezg.download/pd.moex10st.mm.entrezg_25.0.0.tar.gz", repos = NULL, type = "source")
#install.packages("http://mbni.org/customcdf/25.0.0/entrezg.download/mogene10stmmentrezg.db_25.0.0.tar.gz", repos = NULL, type = "source")
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
expb <- exp[,c(1:3, 49:51)]
gse <- getGEO("GSE31997", AnnotGPL = TRUE)
pd <- pData(gse$GSE31997_series_matrix.txt.gz)
NIblb= as.character(pd$title[1:3])
NIC= as.character(pd$title[49:51])
pd.desNIb <- NIblb %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$")%>%
  as.data.frame() %>%
  unite("V2_V3", c("V2","V3"),sep = " " ) %>%
  dplyr::select(V2_V3,V4, V7, V8) %>%
  setNames(c("Genotype","Strain", "Time", "Replicate")) %>% 
  mutate(Genotype = gsub("Non Infected","NIBlb6T0",Genotype),
         Genotype = as.character(Genotype))

pd.desNIC <- NIC %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
  as.data.frame() %>% 
  unite("V2_V3", c("V2","V3"),sep = " " ) %>%
  dplyr::select(V2_V3,V4, V7, V8) %>%
  setNames(c("Genotype", "Strain", "Time", "Replicate")) %>% 
  mutate(Genotype = gsub("Non Infected","NIC57T0",Genotype),
         Genotype = as.character(Genotype))
pd.desb<- bind_rows(pd.desNIb,pd.desNIC)
rownames(pd.desb) <- pd$geo_accession[c(1:3, 49:51)]
Genotype = pd.desb$Genotype
eset_filtereds <- exprs(eset_filtered)[,c(1:3, 49:51)]
require(limma)
designMat<- model.matrix(~0+Genotype)
dim(designMat)
colnames(designMat) <- c("NIBlb6T0","NIC57T0")
cont.matrix <- makeContrasts(NIBlb6vsNIC57T0 = NIBlb6T0-NIC57T0,levels=designMat)
print(cont.matrix)
fit<-lmFit(expb, designMat)
fit.main<-contrasts.fit(fit, cont.matrix)
fit.main<-eBayes(fit.main)
class(fit.main)
topTab_NIBvsNIC <- topTable (fit.main, number=nrow(fit.main), coef="NIBlb6vsNIC57T0", adjust="fdr") 
topTab_NIBvsNIC_filt <- topTab_NIBvsNIC %>% filter(adj.P.Val < 0.05)
NIBvsNIC_filt_all<- topTab_NIBvsNIC_filt %>% filter(-1 > logFC  | logFC > 1)
nrow(NIBvsNIC_filt_all)
NIBvsNIC_filt_up <- topTab_NIBvsNIC_filt %>% filter( logFC > 1)
NIBvsNIC_filt_up
nrow(NIBvsNIC_filt_up)
NIBvsNIC_filt_down <- topTab_NIBvsNIC_filt %>% filter( logFC < -1)
nrow(NIBvsNIC_filt_down)
write.csv(NIBvsNIC_filt_all, "NIBvsNIC_filt_all_T0.csv")
write.csv(NIBvsNIC_filt_up, "NIBvsNIC_filt_up_T0.csv")
write.csv(NIBvsNIC_filt_down, "NIBvsNIC_filt_down_T0.csv")
require(devtools)
source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2entrez.R")
require(biomaRt)
mart.mm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "mmusculus_gene_ensembl",
                   host="useast.ensembl.org")
NIBvsNIC_filt_all_annot   <- probe2entrez(NIBvsNIC_filt_all, mart.mm)
NIBvsNIC_filt_down_annot   <- probe2entrez(NIBvsNIC_filt_down, mart.mm)
NIBvsNIC_filt_up_annot   <- probe2entrez(NIBvsNIC_filt_up, mart.mm)
write.csv(NIBvsNIC_filt_all_annot, "NIBvsNIC_filt_all_annot_T0.csv")
write.csv(NIBvsNIC_filt_up_annot, "NIBvsNIC_filt_up_annot_T0.csv")
write.csv(NIBvsNIC_filt_down_annot, "NIBvsNIC_filt_down_annot_T0.csv")
require(mogene10sttranscriptcluster.db)
geneSymbols <- AnnotationDbi::select(mogene10sttranscriptcluster.db, rownames(topTab_NIBvsNIC_filt), c("SYMBOL"))
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
#### T1 ####
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
#install.packages("http://mbni.org/customcdf/25.0.0/entrezg.download/pd.moex10st.mm.entrezg_25.0.0.tar.gz", repos = NULL, type = "source")
#install.packages("http://mbni.org/customcdf/25.0.0/entrezg.download/mogene10stmmentrezg.db_25.0.0.tar.gz", repos = NULL, type = "source")
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
expb <- exp[,c(4:6, 52:54)]
gse <- getGEO("GSE31997", AnnotGPL = TRUE)
pd <- pData(gse$GSE31997_series_matrix.txt.gz)
NIblb= as.character(pd$title[4:6])
NIC= as.character(pd$title[52:54])
pd.desNIb <- NIblb %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$")%>%
  as.data.frame() %>%
  unite("V2_V3", c("V2","V3"),sep = " " ) %>%
  dplyr::select(V2_V3,V4, V7, V8) %>%
  setNames(c("Genotype","Strain", "Time", "Replicate")) %>% 
  mutate(Genotype = gsub("Non Infected","NIBlb6T1",Genotype),
         Genotype = as.character(Genotype))

pd.desNIC <- NIC %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
  as.data.frame() %>% 
  unite("V2_V3", c("V2","V3"),sep = " " ) %>%
  dplyr::select(V2_V3,V4, V7, V8) %>%
  setNames(c("Genotype", "Strain", "Time", "Replicate")) %>% 
  mutate(Genotype = gsub("Non Infected","NIC57T1",Genotype),
         Genotype = as.character(Genotype))
pd.desb<- bind_rows(pd.desNIb,pd.desNIC)
rownames(pd.desb) <- pd$geo_accession[c(4:6, 52:54)]
Genotype = pd.desb$Genotype
eset_filtereds <- exprs(eset_filtered)[,c(4:6, 52:54)]
require(limma)
designMat<- model.matrix(~0+Genotype)
dim(designMat)
colnames(designMat) <- c("NIBlb6T1","NIC57T1")
cont.matrix <- makeContrasts(NIBlb6vsNIC57T1 = NIBlb6T1-NIC57T1,levels=designMat)
print(cont.matrix)
fit<-lmFit(expb, designMat)
fit.main<-contrasts.fit(fit, cont.matrix)
fit.main<-eBayes(fit.main)
class(fit.main)
topTab_NIBvsNIC <- topTable (fit.main, number=nrow(fit.main), coef="NIBlb6vsNIC57T1", adjust="fdr") 
topTab_NIBvsNIC_filt <- topTab_NIBvsNIC %>% filter(adj.P.Val < 0.05)
NIBvsNIC_filt_up <- topTab_NIBvsNIC_filt %>% filter( logFC > 1)
NIBvsNIC_filt_up
nrow(NIBvsNIC_filt_up)
NIBvsNIC_filt_down <- topTab_NIBvsNIC_filt %>% filter( logFC < -1)
nrow(NIBvsNIC_filt_down)
write.csv(NIBvsNIC_filt_up, "NIBvsNIC_filt_up_T1.csv")
write.csv(NIBvsNIC_filt_down, "NIBvsNIC_filt_down_T1.csv")
require(devtools)
source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2entrez.R")
require(biomaRt)
mart.mm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "mmusculus_gene_ensembl",
                   host="useast.ensembl.org")
NIBvsNIC_filt_down_annot   <- probe2entrez(NIBvsNIC_filt_down, mart.mm)
NIBvsNIC_filt_up_annot   <- probe2entrez(NIBvsNIC_filt_up, mart.mm)
write.csv(NIBvsNIC_filt_up_annot, "NIBvsNIC_filt_up_annot_T1.csv")
write.csv(NIBvsNIC_filt_down_annot, "NIBvsNIC_filt_down_annot_T1.csv")
require(mogene10sttranscriptcluster.db)
geneSymbols <- select(mogene10sttranscriptcluster.db, rownames(topTab_NIBvsNIC_filt), c("SYMBOL"))
SYMBOLS<- geneSymbols$SYMBOL
volcanoplot(fit.main, coef=1, highlight=10, names=SYMBOLS, 
            main=paste("Differentially expressed genes", colnames(cont.matrix)[1], sep="\n"))
abline(v=c(-1,1))

#### T3 ####
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
#install.packages("http://mbni.org/customcdf/25.0.0/entrezg.download/pd.moex10st.mm.entrezg_25.0.0.tar.gz", repos = NULL, type = "source")
#install.packages("http://mbni.org/customcdf/25.0.0/entrezg.download/mogene10stmmentrezg.db_25.0.0.tar.gz", repos = NULL, type = "source")
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
expb <- exp[,c(7:9, 55:57)]
gse <- getGEO("GSE31997", AnnotGPL = TRUE)
pd <- pData(gse$GSE31997_series_matrix.txt.gz)
NIblb= as.character(pd$title[7:9])
NIC= as.character(pd$title[55:57])
pd.desNIb <- NIblb %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$")%>%
  as.data.frame() %>%
  unite("V2_V3", c("V2","V3"),sep = " " ) %>%
  dplyr::select(V2_V3,V4, V7, V8) %>%
  setNames(c("Genotype","Strain", "Time", "Replicate")) %>% 
  mutate(Genotype = gsub("Non Infected","NIBlb6T3",Genotype),
         Genotype = as.character(Genotype))

pd.desNIC <- NIC %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
  as.data.frame() %>% 
  unite("V2_V3", c("V2","V3"),sep = " " ) %>%
  dplyr::select(V2_V3,V4, V7, V8) %>%
  setNames(c("Genotype", "Strain", "Time", "Replicate")) %>% 
  mutate(Genotype = gsub("Non Infected","NIC57T3",Genotype),
         Genotype = as.character(Genotype))
pd.desb<- bind_rows(pd.desNIb,pd.desNIC)
rownames(pd.desb) <- pd$geo_accession[c(7:9, 55:57)]
Genotype = pd.desb$Genotype
eset_filtereds <- exprs(eset_filtered)[,c(7:9, 55:57)]
require(limma)
designMat<- model.matrix(~0+Genotype)
dim(designMat)
colnames(designMat) <- c("NIBlb6T3","NIC57T3")
cont.matrix <- makeContrasts(NIBlb6vsNIC57T3 = NIBlb6T3-NIC57T3,levels=designMat)
print(cont.matrix)
fit<-lmFit(expb, designMat)
fit.main<-contrasts.fit(fit, cont.matrix)
fit.main<-eBayes(fit.main)
class(fit.main)
topTab_NIBvsNIC <- topTable (fit.main, number=nrow(fit.main), coef="NIBlb6vsNIC57T3", adjust="fdr") 
topTab_NIBvsNIC_filt <- topTab_NIBvsNIC %>% filter(adj.P.Val < 0.05)
NIBvsNIC_filt_up <- topTab_NIBvsNIC_filt %>% filter( logFC > 1)
NIBvsNIC_filt_up
nrow(NIBvsNIC_filt_up)
NIBvsNIC_filt_down <- topTab_NIBvsNIC_filt %>% filter( logFC < -1)
nrow(NIBvsNIC_filt_down)
write.csv(NIBvsNIC_filt_up, "NIBvsNIC_filt_up_T3.csv")
write.csv(NIBvsNIC_filt_down, "NIBvsNIC_filt_down_T3.csv")
require(devtools)
source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2entrez.R")
require(biomaRt)
mart.mm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "mmusculus_gene_ensembl",
                   host="useast.ensembl.org")
NIBvsNIC_filt_down_annot   <- probe2entrez(NIBvsNIC_filt_down, mart.mm)
NIBvsNIC_filt_up_annot   <- probe2entrez(NIBvsNIC_filt_up, mart.mm)
write.csv(NIBvsNIC_filt_up_annot, "NIBvsNIC_filt_up_annot_T3.csv")
write.csv(NIBvsNIC_filt_down_annot, "NIBvsNIC_filt_down_annot_T3.csv")
require(mogene10sttranscriptcluster.db)
geneSymbols <- select(mogene10sttranscriptcluster.db, rownames(topTab_NIBvsNIC_filt), c("SYMBOL"))
SYMBOLS<- geneSymbols$SYMBOL
volcanoplot(fit.main, coef=1, highlight=10, names=SYMBOLS, 
            main=paste("Differentially expressed genes", colnames(cont.matrix)[1], sep="\n"))
abline(v=c(-1,1))

#### T6 ####
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
#install.packages("http://mbni.org/customcdf/25.0.0/entrezg.download/pd.moex10st.mm.entrezg_25.0.0.tar.gz", repos = NULL, type = "source")
#install.packages("http://mbni.org/customcdf/25.0.0/entrezg.download/mogene10stmmentrezg.db_25.0.0.tar.gz", repos = NULL, type = "source")
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
expb <- exp[,c(10:12, 58:60)]
gse <- getGEO("GSE31997", AnnotGPL = TRUE)
pd <- pData(gse$GSE31997_series_matrix.txt.gz)
NIblb= as.character(pd$title[10:12])
NIC= as.character(pd$title[58:60])
pd.desNIb <- NIblb %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$")%>%
  as.data.frame() %>%
  unite("V2_V3", c("V2","V3"),sep = " " ) %>%
  dplyr::select(V2_V3,V4, V7, V8) %>%
  setNames(c("Genotype","Strain", "Time", "Replicate")) %>% 
  mutate(Genotype = gsub("Non Infected","NIBlb6T1",Genotype),
         Genotype = as.character(Genotype))

pd.desNIC <- NIC %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
  as.data.frame() %>% 
  unite("V2_V3", c("V2","V3"),sep = " " ) %>%
  dplyr::select(V2_V3,V4, V7, V8) %>%
  setNames(c("Genotype", "Strain", "Time", "Replicate")) %>% 
  mutate(Genotype = gsub("Non Infected","NIC57T1",Genotype),
         Genotype = as.character(Genotype))
pd.desb<- bind_rows(pd.desNIb,pd.desNIC)
rownames(pd.desb) <- pd$geo_accession[c(10:12, 58:60)]
Genotype = pd.desb$Genotype
eset_filtereds <- exprs(eset_filtered)[,c(10:12, 58:60)]
require(limma)
designMat<- model.matrix(~0+Genotype)
dim(designMat)
colnames(designMat) <- c("NIBlb6T6","NIC57T6")
cont.matrix <- makeContrasts(NIBlb6vsNIC57T6 = NIBlb6T6-NIC57T6,levels=designMat)
print(cont.matrix)
fit<-lmFit(expb, designMat)
fit.main<-contrasts.fit(fit, cont.matrix)
fit.main<-eBayes(fit.main)
class(fit.main)
topTab_NIBvsNIC <- topTable (fit.main, number=nrow(fit.main), coef="NIBlb6vsNIC57T6", adjust="fdr") 
topTab_NIBvsNIC_filt <- topTab_NIBvsNIC %>% filter(adj.P.Val < 0.05)
NIBvsNIC_filt_up <- topTab_NIBvsNIC_filt %>% filter( logFC > 1)
NIBvsNIC_filt_up
nrow(NIBvsNIC_filt_up)
NIBvsNIC_filt_down <- topTab_NIBvsNIC_filt %>% filter( logFC < -1)
nrow(NIBvsNIC_filt_down)
write.csv(NIBvsNIC_filt_up, "NIBvsNIC_filt_up_T6.csv")
write.csv(NIBvsNIC_filt_down, "NIBvsNIC_filt_down_T6.csv")
require(devtools)
source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2entrez.R")
require(biomaRt)
mart.mm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "mmusculus_gene_ensembl",
                   host="useast.ensembl.org")
NIBvsNIC_filt_down_annot   <- probe2entrez(NIBvsNIC_filt_down, mart.mm)
NIBvsNIC_filt_up_annot   <- probe2entrez(NIBvsNIC_filt_up, mart.mm)
write.csv(NIBvsNIC_filt_up_annot, "NIBvsNIC_filt_up_annot_T6.csv")
write.csv(NIBvsNIC_filt_down_annot, "NIBvsNIC_filt_down_annot_T6.csv")
require(mogene10sttranscriptcluster.db)
geneSymbols <- select(mogene10sttranscriptcluster.db, rownames(topTab_NIBvsNIC_filt), c("SYMBOL"))
SYMBOLS<- geneSymbols$SYMBOL
volcanoplot(fit.main, coef=1, highlight=10, names=SYMBOLS, 
            main=paste("Differentially expressed genes", colnames(cont.matrix)[1], sep="\n"))
abline(v=c(-1,1))

#### T12 ####
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
#install.packages("http://mbni.org/customcdf/25.0.0/entrezg.download/pd.moex10st.mm.entrezg_25.0.0.tar.gz", repos = NULL, type = "source")
#install.packages("http://mbni.org/customcdf/25.0.0/entrezg.download/mogene10stmmentrezg.db_25.0.0.tar.gz", repos = NULL, type = "source")
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
expb <- exp[,c(13:15, 61:63)]
gse <- getGEO("GSE31997", AnnotGPL = TRUE)
pd <- pData(gse$GSE31997_series_matrix.txt.gz)
NIblb= as.character(pd$title[13:15])
NIC= as.character(pd$title[61:63])
pd.desNIb <- NIblb %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$")%>%
  as.data.frame() %>%
  unite("V2_V3", c("V2","V3"),sep = " " ) %>%
  dplyr::select(V2_V3,V4, V7, V8) %>%
  setNames(c("Genotype","Strain", "Time", "Replicate")) %>% 
  mutate(Genotype = gsub("Non Infected","NIBlb6T12",Genotype),
         Genotype = as.character(Genotype))

pd.desNIC <- NIC %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
  as.data.frame() %>% 
  unite("V2_V3", c("V2","V3"),sep = " " ) %>%
  dplyr::select(V2_V3,V4, V7, V8) %>%
  setNames(c("Genotype", "Strain", "Time", "Replicate")) %>% 
  mutate(Genotype = gsub("Non Infected","NIC57T12",Genotype),
         Genotype = as.character(Genotype))
pd.desb<- bind_rows(pd.desNIb,pd.desNIC)
rownames(pd.desb) <- pd$geo_accession[c(13:15, 61:63)]
Genotype = pd.desb$Genotype
eset_filtereds <- exprs(eset_filtered)[,c(13:15, 61:63)]
require(limma)
designMat<- model.matrix(~0+Genotype)
dim(designMat)
colnames(designMat) <- c("NIBlb6T12","NIC57T12")
cont.matrix <- makeContrasts(NIBlb6vsNIC57T12 = NIBlb6T12-NIC57T12,levels=designMat)
print(cont.matrix)
fit<-lmFit(expb, designMat)
fit.main<-contrasts.fit(fit, cont.matrix)
fit.main<-eBayes(fit.main)
class(fit.main)
topTab_NIBvsNIC <- topTable (fit.main, number=nrow(fit.main), coef="NIBlb6vsNIC57T12", adjust="fdr") 
topTab_NIBvsNIC_filt <- topTab_NIBvsNIC %>% filter(adj.P.Val < 0.05)
NIBvsNIC_filt_up <- topTab_NIBvsNIC_filt %>% filter( logFC > 1)
NIBvsNIC_filt_up
nrow(NIBvsNIC_filt_up)
NIBvsNIC_filt_down <- topTab_NIBvsNIC_filt %>% filter( logFC < -1)
nrow(NIBvsNIC_filt_down)
write.csv(NIBvsNIC_filt_up, "NIBvsNIC_filt_up_T12.csv")
write.csv(NIBvsNIC_filt_down, "NIBvsNIC_filt_down_T12.csv")
require(devtools)
source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2entrez.R")
require(biomaRt)
mart.mm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "mmusculus_gene_ensembl",
                   host="useast.ensembl.org")
NIBvsNIC_filt_down_annot   <- probe2entrez(NIBvsNIC_filt_down, mart.mm)
NIBvsNIC_filt_up_annot   <- probe2entrez(NIBvsNIC_filt_up, mart.mm)
write.csv(NIBvsNIC_filt_up_annot, "NIBvsNIC_filt_up_annot_T12.csv")
write.csv(NIBvsNIC_filt_down_annot, "NIBvsNIC_filt_down_annot_T12.csv")
require(mogene10sttranscriptcluster.db)
geneSymbols <- select(mogene10sttranscriptcluster.db, rownames(topTab_NIBvsNIC_filt), c("SYMBOL"))
SYMBOLS<- geneSymbols$SYMBOL
volcanoplot(fit.main, coef=1, highlight=10, names=SYMBOLS, 
            main=paste("Differentially expressed genes", colnames(cont.matrix)[1], sep="\n"))
abline(v=c(-1,1))

#### T24 ####
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
#install.packages("http://mbni.org/customcdf/25.0.0/entrezg.download/pd.moex10st.mm.entrezg_25.0.0.tar.gz", repos = NULL, type = "source")
#install.packages("http://mbni.org/customcdf/25.0.0/entrezg.download/mogene10stmmentrezg.db_25.0.0.tar.gz", repos = NULL, type = "source")
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
expb <- exp[,c(16:18, 64:66)]
gse <- getGEO("GSE31997", AnnotGPL = TRUE)
pd <- pData(gse$GSE31997_series_matrix.txt.gz)
NIblb= as.character(pd$title[16:18])
NIC= as.character(pd$title[64:66])
pd.desNIb <- NIblb %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$")%>%
  as.data.frame() %>%
  unite("V2_V3", c("V2","V3"),sep = " " ) %>%
  dplyr::select(V2_V3,V4, V7, V8) %>%
  setNames(c("Genotype","Strain", "Time", "Replicate")) %>% 
  mutate(Genotype = gsub("Non Infected","NIBlb6T24",Genotype),
         Genotype = as.character(Genotype))

pd.desNIC <- NIC %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
  as.data.frame() %>% 
  unite("V2_V3", c("V2","V3"),sep = " " ) %>%
  dplyr::select(V2_V3,V4, V7, V8) %>%
  setNames(c("Genotype", "Strain", "Time", "Replicate")) %>% 
  mutate(Genotype = gsub("Non Infected","NIC57T24",Genotype),
         Genotype = as.character(Genotype))
pd.desb<- bind_rows(pd.desNIb,pd.desNIC)
rownames(pd.desb) <- pd$geo_accession[c(16:18, 64:66)]
Genotype = pd.desb$Genotype
eset_filtereds <- exprs(eset_filtered)[,c(16:18, 64:66)]
require(limma)
designMat<- model.matrix(~0+Genotype)
dim(designMat)
colnames(designMat) <- c("NIBlb6T24","NIC57T24")
cont.matrix <- makeContrasts(NIBlb6vsNIC57T24 = NIBlb6T24-NIC57T24,levels=designMat)
print(cont.matrix)
fit<-lmFit(expb, designMat)
fit.main<-contrasts.fit(fit, cont.matrix)
fit.main<-eBayes(fit.main)
class(fit.main)
topTab_NIBvsNIC <- topTable (fit.main, number=nrow(fit.main), coef="NIBlb6vsNIC57T24", adjust="fdr") 
topTab_NIBvsNIC_filt <- topTab_NIBvsNIC %>% filter(adj.P.Val < 0.05)
NIBvsNIC_filt_all<- topTab_NIBvsNIC_filt %>% filter(-1 > logFC  & logFC < 1)
NIBvsNIC_filt_up <- topTab_NIBvsNIC_filt %>% filter( logFC > 1)
NIBvsNIC_filt_up
nrow(NIBvsNIC_filt_up)
NIBvsNIC_filt_down <- topTab_NIBvsNIC_filt %>% filter( logFC < -1)
nrow(NIBvsNIC_filt_down)
write.csv(NIBvsNIC_filt_up, "NIBvsNIC_filt_up_T24.csv")
write.csv(NIBvsNIC_filt_down, "NIBvsNIC_filt_down_T24.csv")
require(devtools)
source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2entrez.R")
require(biomaRt)
mart.mm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "mmusculus_gene_ensembl",
                   host="useast.ensembl.org")
NIBvsNIC_filt_all_annot   <- probe2entrez(NIBvsNIC_filt_all, mart.mm)
NIBvsNIC_filt_down_annot   <- probe2entrez(NIBvsNIC_filt_down, mart.mm)
NIBvsNIC_filt_up_annot   <- probe2entrez(NIBvsNIC_filt_up, mart.mm)
write.csv(NIBvsNIC_filt_up_annot, "NIBvsNIC_filt_up_annot_T24.csv")
write.csv(NIBvsNIC_filt_down_annot, "NIBvsNIC_filt_down_annot_T24.csv")
require(mogene10sttranscriptcluster.db)
geneSymbols <- select(mogene10sttranscriptcluster.db, rownames(topTab_NIBvsNIC_filt), c("SYMBOL"))
SYMBOLS<- geneSymbols$SYMBOL
volcanoplot(fit.main, coef=1, highlight=10, names=SYMBOLS, 
            main=paste("Differentially expressed genes", colnames(cont.matrix)[1], sep="\n"))
abline(v=c(-1,1))
#################### Enrichment ##########################
library(goseq)
supportedOrganisms() %>% filter(str_detect(Genome, "mm"))


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
logFC <- NIBvsNIC_filt_down_annot$logFC
names(logFC) <- NIBvsNIC_filt_all_annot$entrez_id
pathview(gene.data = logFC, 
         pathway.id = "mmu04610", 
         species = "mmu", 
         limit = list(gene=5, cpd=1))
