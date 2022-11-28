# Author: Hédia Tnani
# Date: 04/03/2020
# Usage:MasigLeish(mouse = "Balb", comp=c("IvsC","PvsC","KPvsC","PvsKP")) or MasigLeish(mouse = "c57", comp=c("IvsC","PvsC","KPvsC","PvsKP"))
#################################################################################################################################
MasigLeish = function (mouse, comp){
  options(warn = -1)
  cat ("Install required packages")
  packages <- c("tidyverse", "doParallel","foreign","itertools", "devtools", "dplyr")
  for (pkg in packages) {
    if (!pkg %in% rownames(installed.packages())) {
      message("Installing missing package: ", pkg)
      install.packages(pkg)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  requiredPackages=c("oligo","GEOquery","oligo","Biobase","affy","mogene10sttranscriptcluster.db", "maSigPro", "biomaRt","qvalue","RCurl")
  for (pkg_name in requiredPackages) {
    if (!pkg_name %in% rownames(installed.packages())) {
      message("Installing missing Bioconductor package: ", pkg_name)
      if (!requireNamespace("BiocManager", quietly = TRUE)){
        install.packages("BiocManager")
        BiocManager::install(version = "3.10")
        BiocManager::install(pkg_name)
      }}
    suppressPackageStartupMessages(library(pkg_name, character.only = TRUE))
  }
  #### Balb ####
  if (mouse== "Balb"){# Balb
    message ("loading packages ....")
    require(GEOquery)
    require(oligo)
    require(mogene10sttranscriptcluster.db)
    require(RCurl)
    require(foreign)
    require(tidyverse)
    setwd(getwd())
    if (comp== "IvsC"){
      if (!dir.exists(paste0(getwd(),"/gse31997data"))){
        message("downloading data ...")
        getGEOSuppFiles("GSE31997")
        untar("GSE31997/GSE31997_RAW.tar", exdir="./gse31997data")
        cels <- list.files("gse31997data/", pattern = "[gz]")
        message("unzipping files ...")
        sapply(paste("gse31997data", cels, sep="/"), gunzip)
        celFiles <- list.celfiles("./gse31997data",full.names = TRUE)
        message("getting targets data from github ...")
        url <- "https://raw.githubusercontent.com/HediaTnani/Leishmania/master/targets.csv"
        targets <- getURL(url)                
        my.targets <- read.AnnotatedDataFrame(textConnection(targets),header = TRUE, row.names = 1, sep=",")
        message("Listing celfiles...")
        rawData <- read.celfiles(celFiles, phenoData = my.targets)
        colnames(rawData) <-rownames(pData(rawData))<-my.targets@data$ShortName
        eset_rma <- oligo::rma(rawData) # specify rma from package oligo or affy
        exp <- data.frame(exprs(eset_rma))
        expb <- exp[,4:48]
        message("Viewing expression data ...")
        expb %>%
          as.data.frame() %>%
          head()%>%
          dplyr::select(1:4)
        message("Getting phenotypic data ...")
        gse <- getGEO("GSE31997", AnnotGPL = TRUE)
        pd <- pData(gse$GSE31997_series_matrix.txt.gz)
        NIb= as.character(pd$title[4:18])
        Ib= as.character(pd$title[19:48])
        message("Printing phenotypic data ...")
        print(NIb)
        print(Ib)
        message("Creating a design matrix ...")
        pd.desNIb <- NIb %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>%
          as.data.frame() %>%
          unite("V2_V3", c("V2","V3"),sep = " " ) %>%
          dplyr::select(V2_V3,V7, V8) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate(Control = ifelse(Genotype == "Non Infected", 1, 0),
                 Genotype = gsub("[[:punct:]]*", "", Genotype),
                 Genotype = as.character(Genotype),
                 Infected = ifelse(Genotype == "P"|Genotype == "KP", 1, 0),
                 Time = gsub("T", "", Time),
                 Time = as.numeric(Time))
        pd.desIb <- Ib %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>%
          as.data.frame() %>%
          dplyr::select(V8,V10, V11) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate( Control = ifelse(Genotype == "Non Infected", 1, 0),
                  Genotype = gsub("[[:punct:]]*", "", Genotype),
                  Genotype = as.character(Genotype),
                  Infected = ifelse(Genotype == "P"|Genotype == "KP", 1, 0),
                  Time = gsub("T", "", Time),
                  Time = as.numeric(Time))
        pd.desb<- bind_rows(pd.desNIb,pd.desIb)
        pd.desb$Replicate<- rep(1:15, 1, each = 3)
        DT::datatable(head(pd.desb), caption = "Design matrix")
        rownames(pd.desb) <- pd$geo_accession[4:48]
        require(maSigPro)
        designb <- make.design.matrix(pd.desb[, c(2,3,4,5)], degree = 2)
        designb$edesign %>% as.data.frame() %>% DT::datatable()
        message("Defining the regression model...")
        designb$groups.vector
        message("sourcing functions from github...")
        require(doParallel)
        require(itertools)
        require(devtools)
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/Tfitp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/pvectorp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/position.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2sym.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2entrez.R")
        message("Fitting the regression model...")
        fitb<- p.vectorp(expb, designb,Q = 0.1, MT.adjust = "BH", min.obs = 20)
        tstepb <- T.fitp(fitb, step.method = "backward", alfa = 0.05)
        sigsb <- get.siggenes(tstepb, rsq = 0.6, vars = "groups")
        sigsb %>% glimpse()
        names(sigsb)
        names(sigsb$sig.genes)
        suma2Venn(sigsb$summary[,c(1:2)])
        clusters=see.genes(sigsb$sig.genes$InfectedvsControl, show.fit = T, dis=designb$dis, cluster.method = "hclust", cluster.data = 1, k.mclust = TRUE)
        df=clusters$cut %>% as.data.frame() 
        write.csv(df, "IvsC_B_Geneclusters.csv")
        IvsC   <- sigsb$sig.genes$InfectedvsControl$sig.pvalues
        message ("from probes to entrez id ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "mmusculus_gene_ensembl",
                           host = 'www.ensembl.org')
        IvsC   <- probe2entrez(IvsC, mart.mm)
        #IvsC1 <- data.frame(IvsC, row.names = 1)
        write.csv(IvsC, "IvsC_B_v0.csv")
        message ("from probes to genes symbols ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ensembl",
                         dataset = "mmusculus_gene_ensembl",
                         host="https://dec2021.archive.ensembl.org")
        mart.hs = useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                        host="https://dec2021.archive.ensembl.org")
        IvsC   <- probe2sym(IvsC, mart.mm,mart.hs)
        message ("writing csv files ....")
        write.csv(IvsC, "IvsC_B.csv",row.names=TRUE)
        message ("Analysis IvsC for Balb completed successfully ....")
      }else{
        celFiles <- list.celfiles("./gse31997data",full.names = TRUE)
        message("getting targets data from github ...")
        url <- "https://raw.githubusercontent.com/HediaTnani/Leishmania/master/targets.csv"
        targets <- getURL(url)                
        my.targets <- read.AnnotatedDataFrame(textConnection(targets),header = TRUE, row.names = 1, sep=",")
        message("Listing celfiles...")
        rawData <- read.celfiles(celFiles, phenoData = my.targets)
        colnames(rawData) <-rownames(pData(rawData))<-my.targets@data$ShortName
        eset_rma <- oligo::rma(rawData) # specify rma from package oligo or affy
        exp <- data.frame(exprs(eset_rma))
        expb <- exp[,4:48]
        message("Viewing expression data ...")
        expb %>%
          as.data.frame() %>%
          head()%>%
          dplyr::select(1:4)
        message("Getting phenotypic data ...")
        gse <- getGEO("GSE31997", AnnotGPL = TRUE)
        pd <- pData(gse$GSE31997_series_matrix.txt.gz)
        NIb= as.character(pd$title[4:18])
        Ib= as.character(pd$title[19:48])
        message("Printing phenotypic data ...")
        print(NIb)
        print(Ib)
        message("Creating a design matrix ...")
        pd.desNIb <- NIb %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>%
          as.data.frame() %>%
          unite("V2_V3", c("V2","V3"),sep = " " ) %>%
          dplyr::select(V2_V3,V7, V8) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate(Control = ifelse(Genotype == "Non Infected", 1, 0),
                 Genotype = gsub("[[:punct:]]*", "", Genotype),
                 Genotype = as.character(Genotype),
                 Infected = ifelse(Genotype == "P"|Genotype == "KP", 1, 0),
                 Time = gsub("T", "", Time),
                 Time = as.numeric(Time))
        pd.desIb <- Ib %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>%
          as.data.frame() %>%
          dplyr::select(V8,V10, V11) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate( Control = ifelse(Genotype == "Non Infected", 1, 0),
                  Genotype = gsub("[[:punct:]]*", "", Genotype),
                  Genotype = as.character(Genotype),
                  Infected = ifelse(Genotype == "P"|Genotype == "KP", 1, 0),
                  Time = gsub("T", "", Time),
                  Time = as.numeric(Time))
        pd.desb<- bind_rows(pd.desNIb,pd.desIb)
        pd.desb$Replicate<- rep(1:15, 1, each = 3)
        DT::datatable(head(pd.desb), caption = "Design matrix")
        rownames(pd.desb) <- pd$geo_accession[4:48]
        require(maSigPro)
        designb <- make.design.matrix(pd.desb[, c(2,3,4,5)], degree = 2)
        designb$edesign %>% as.data.frame() %>% DT::datatable()
        message("Defining the regression model...")
        designb$groups.vector
        message("sourcing functions from github...")
        require(doParallel)
        require(itertools)
        require(devtools)
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/Tfitp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/pvectorp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/position.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2sym.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2entrez.R")
        message("Fitting the regression model...")
        fitb<- p.vectorp(expb, designb,Q = 0.1, MT.adjust = "BH", min.obs = 20)
        tstepb <- T.fitp(fitb, step.method = "backward", alfa = 0.05)
        sigsb <- get.siggenes(tstepb, rsq = 0.6, vars = "groups")
        sigsb %>% glimpse()
        names(sigsb)
        names(sigsb$sig.genes)
        suma2Venn(sigsb$summary[,c(1:2)])
        clusters=see.genes(sigsb$sig.genes$InfectedvsControl, show.fit = T, dis=designb$dis, cluster.method = "hclust", cluster.data = 1, k.mclust = TRUE)
        df=clusters$cut %>% as.data.frame() 
        write.csv(df, "IvsC_B_Geneclusters.csv")
        IvsC   <- sigsb$sig.genes$InfectedvsControl$sig.pvalues
        message ("from probes to entrez id ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "mmusculus_gene_ensembl",
                           host = 'www.ensembl.org')
        IvsC   <- probe2entrez(IvsC, mart.mm)
        write.csv(IvsC, "IvsC_B_v0.csv")
        message ("from probes to genes symbols ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ensembl",
                         dataset = "mmusculus_gene_ensembl",
                         host="https://dec2021.archive.ensembl.org")
        mart.hs = useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                        host="https://dec2021.archive.ensembl.org")
        IvsC   <- probe2sym(IvsC, mart.mm,mart.hs)
        #IvsC1 <- data.frame(IvsC, row.names = 1)
        message ("writing csv files ....")
        write.csv(IvsC, "IvsC_B.csv",row.names=F)
        message ("Analysis IvsC for Balb completed successfully ....")
      }#else
    }else if (comp== "PvsC"){
      if (!dir.exists(paste0(getwd(),"/gse31997data"))){
        message("downloading data ...")
        getGEOSuppFiles("GSE31997")
        untar("GSE31997/GSE31997_RAW.tar", exdir="./gse31997data")
        cels <- list.files("gse31997data/", pattern = "[gz]")
        message("unzipping files ...")
        sapply(paste("gse31997data", cels, sep="/"), gunzip)
        celFiles <- list.celfiles("./gse31997data",full.names = TRUE)
        message("getting targets data from github ...")
        url <- "https://raw.githubusercontent.com/HediaTnani/Leishmania/master/targets.csv"
        targets <- getURL(url)                
        my.targets <- read.AnnotatedDataFrame(textConnection(targets),header = TRUE, row.names = 1, sep=",")
        message("Listing celfiles...")
        rawData <- read.celfiles(celFiles, phenoData = my.targets)
        colnames(rawData) <-rownames(pData(rawData))<-my.targets@data$ShortName
        eset_rma <- oligo::rma(rawData) # specify rma from package oligo or affy
        exp <- data.frame(exprs(eset_rma))
        expb <- exp[,4:33]
        message("Viewing expression data ...")
        expb %>%
          as.data.frame() %>%
          head()%>%
          dplyr::select(1:4)
        message("Getting phenotypic data ...")
        gse <- getGEO("GSE31997", AnnotGPL = TRUE)
        pd <- pData(gse$GSE31997_series_matrix.txt.gz)
        NIb= as.character(pd$title[4:18])
        Ib= as.character(pd$title[19:33])
        message("Printing phenotypic data ...")
        print(NIb)
        print(Ib)
        message("Creating a design matrix ...")
        pd.desNIb <- NIb %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>%
          as.data.frame() %>%
          unite("V2_V3", c("V2","V3"),sep = " " ) %>%
          dplyr::select(V2_V3,V7, V8) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate(Control = ifelse(Genotype == "Non Infected", 1, 0),
                 Genotype = gsub("[[:punct:]]*", "", Genotype),
                 Genotype = as.character(Genotype),
                 P = ifelse(Genotype == "P", 1, 0),
                 Time = gsub("T", "", Time),
                 Time = as.numeric(Time))
        pd.desIb <- Ib %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>%
          as.data.frame() %>%
          dplyr::select(V8,V10, V11) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate( Control = ifelse(Genotype == "Non Infected", 1, 0),
                  Genotype = gsub("[[:punct:]]*", "", Genotype),
                  Genotype = as.character(Genotype),
                  P = ifelse(Genotype == "P", 1, 0),
                  Time = gsub("T", "", Time),
                  Time = as.numeric(Time))
        pd.desb<- bind_rows(pd.desNIb,pd.desIb)
        pd.desb$Replicate<- rep(1:10, 1, each = 3)
        DT::datatable(head(pd.desb), caption = "Design matrix")
        rownames(pd.desb) <- pd$geo_accession[4:33]
        require(maSigPro)
        designb <- make.design.matrix(pd.desb[, c(2,3,4,5)], degree = 2)
        designb$edesign %>% as.data.frame() %>% DT::datatable()
        message("Defining the regression model...")
        designb$groups.vector
        message("sourcing functions from github...")
        require(doParallel)
        require(itertools)
        require(devtools)
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/Tfitp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/pvectorp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/position.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2sym.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2entrez.R")
        message("Fitting the regression model...")
        fitb<- p.vectorp(expb, designb,Q = 0.1, MT.adjust = "BH", min.obs = 20)
        tstepb <- T.fitp(fitb, step.method = "backward", alfa = 0.05)
        sigsb <- get.siggenes(tstepb, rsq = 0.6, vars = "groups")
        sigsb %>% glimpse()
        names(sigsb)
        names(sigsb$sig.genes)
        suma2Venn(sigsb$summary[,c(1:2)])
        clusters=see.genes(sigsb$sig.genes$PvsControl, show.fit = T, dis=designb$dis, cluster.method = "hclust", cluster.data = 1, k.mclust = TRUE)
        df=clusters$cut %>% as.data.frame() 
        write.csv(df, "PvsC_B_Geneclusters.csv")
        PvsC   <- sigsb$sig.genes$PvsControl$sig.pvalues
        message ("from probes to entrez id ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "mmusculus_gene_ensembl",
                           host = 'www.ensembl.org')
        PvsC    <- probe2entrez(PvsC , mart.mm)
        write.csv(PvsC, "PvsC_B_v0.csv")
        message ("from probes to genes symbols ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ensembl",
                         dataset = "mmusculus_gene_ensembl",
                         host="https://dec2021.archive.ensembl.org")
        mart.hs = useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                        host="https://dec2021.archive.ensembl.org")
        PvsC   <- probe2sym(PvsC, mart.mm,mart.hs)
        message ("writing csv files ....")
        write.csv(PvsC, "PvsC_B.csv")
        message ("Analysis PvsC for Balb completed successfully ....")
      }else{
        celFiles <- list.celfiles("./gse31997data",full.names = TRUE)
        message("getting targets data from github ...")
        url <- "https://raw.githubusercontent.com/HediaTnani/Leishmania/master/targets.csv"
        targets <- getURL(url)                
        my.targets <- read.AnnotatedDataFrame(textConnection(targets),header = TRUE, row.names = 1, sep=",")
        message("Listing celfiles...")
        rawData <- read.celfiles(celFiles, phenoData = my.targets)
        colnames(rawData) <-rownames(pData(rawData))<-my.targets@data$ShortName
        eset_rma <- oligo::rma(rawData) # specify rma from package oligo or affy
        exp <- data.frame(exprs(eset_rma))
        expb <- exp[,4:33]
        message("Viewing expression data ...")
        expb %>%
          as.data.frame() %>%
          head()%>%
          dplyr::select(1:4)
        message("Getting phenotypic data ...")
        gse <- getGEO("GSE31997", AnnotGPL = TRUE)
        pd <- pData(gse$GSE31997_series_matrix.txt.gz)
        NIb= as.character(pd$title[4:18])
        Ib= as.character(pd$title[19:33])
        message("Printing phenotypic data ...")
        print(NIb)
        print(Ib)
        message("Creating a design matrix ...")
        pd.desNIb <- NIb %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>%
          as.data.frame() %>%
          unite("V2_V3", c("V2","V3"),sep = " " ) %>%
          dplyr::select(V2_V3,V7, V8) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate(Control = ifelse(Genotype == "Non Infected", 1, 0),
                 Genotype = gsub("[[:punct:]]*", "", Genotype),
                 Genotype = as.character(Genotype),
                 P = ifelse(Genotype == "P", 1, 0),
                 Time = gsub("T", "", Time),
                 Time = as.numeric(Time))
        pd.desIb <- Ib %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>%
          as.data.frame() %>%
          dplyr::select(V8,V10, V11) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate( Control = ifelse(Genotype == "Non Infected", 1, 0),
                  Genotype = gsub("[[:punct:]]*", "", Genotype),
                  Genotype = as.character(Genotype),
                  P = ifelse(Genotype == "P", 1, 0),
                  Time = gsub("T", "", Time),
                  Time = as.numeric(Time))
        pd.desb<- bind_rows(pd.desNIb,pd.desIb)
        pd.desb$Replicate<- rep(1:10, 1, each = 3)
        DT::datatable(head(pd.desb), caption = "Design matrix")
        rownames(pd.desb) <- pd$geo_accession[4:33]
        require(maSigPro)
        designb <- make.design.matrix(pd.desb[, c(2,3,4,5)], degree = 2)
        designb$edesign %>% as.data.frame() %>% DT::datatable()
        message("Defining the regression model...")
        designb$groups.vector
        message("sourcing functions from github...")
        require(doParallel)
        require(itertools)
        require(devtools)
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/Tfitp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/pvectorp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/position.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2sym.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2entrez.R")
        message("Fitting the regression model...")
        fitb<- p.vectorp(expb, designb,Q = 0.1, MT.adjust = "BH", min.obs = 20)
        tstepb <- T.fitp(fitb, step.method = "backward", alfa = 0.05)
        sigsb <- get.siggenes(tstepb, rsq = 0.6, vars = "groups")
        sigsb %>% glimpse()
        names(sigsb)
        names(sigsb$sig.genes)
        suma2Venn(sigsb$summary[,c(1:2)])
        clusters=see.genes(sigsb$sig.genes$PvsControl, show.fit = T, dis=designb$dis, cluster.method = "hclust", cluster.data = 1, k.mclust = TRUE)
        df=clusters$cut %>% as.data.frame() 
        write.csv(df, "PvsC_B_Geneclusters.csv")
        PvsC   <- sigsb$sig.genes$PvsControl$sig.pvalues
        message ("from probes to entrez id ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "mmusculus_gene_ensembl",
                           host = 'www.ensembl.org')
        PvsC   <- probe2entrez(PvsC, mart.mm)
        write.csv(PvsC, "PvsC_B_v0.csv")
        message ("from probes to genes symbols ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ensembl",
                         dataset = "mmusculus_gene_ensembl",
                         host="https://dec2021.archive.ensembl.org")
        mart.hs = useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                        host="https://dec2021.archive.ensembl.org")
        PvsC   <- probe2sym(PvsC, mart.mm, mart.hs)
        message ("writing csv files ....")
        write.csv(PvsC, "PvsC_B.csv")
        message ("Analysis PvsC for Balb completed successfully ....")
      }#else
    }else if (comp== "KPvsC"){
      if (!dir.exists(paste0(getwd(),"/gse31997data"))){
        message("downloading data ...")
        getGEOSuppFiles("GSE31997")
        untar("GSE31997/GSE31997_RAW.tar", exdir="./gse31997data")
        cels <- list.files("gse31997data/", pattern = "[gz]")
        message("unzipping files ...")
        sapply(paste("gse31997data", cels, sep="/"), gunzip)
        celFiles <- list.celfiles("./gse31997data",full.names = TRUE)
        message("getting targets data from github ...")
        url <- "https://raw.githubusercontent.com/HediaTnani/Leishmania/master/targets.csv"
        targets <- getURL(url)                
        my.targets <- read.AnnotatedDataFrame(textConnection(targets),header = TRUE, row.names = 1, sep=",")
        message("Listing celfiles...")
        rawData <- read.celfiles(celFiles, phenoData = my.targets)
        colnames(rawData) <-rownames(pData(rawData))<-my.targets@data$ShortName
        eset_rma <- oligo::rma(rawData) # specify rma from package oligo or affy
        exp <- data.frame(exprs(eset_rma))
        expb <- exp[,c(4:18,34:48)]
        message("Viewing expression data ...")
        expb %>%
          as.data.frame() %>%
          head()%>%
          dplyr::select(1:4)
        message("Getting phenotypic data ...")
        gse <- getGEO("GSE31997", AnnotGPL = TRUE)
        pd <- pData(gse$GSE31997_series_matrix.txt.gz)
        NIb= as.character(pd$title[4:18])
        Ib= as.character(pd$title[34:48])
        message("Printing phenotypic data ...")
        print(NIb)
        print(Ib)
        message("Creating a design matrix ...")
        pd.desNIb <- NIb %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>%
          as.data.frame() %>%
          unite("V2_V3", c("V2","V3"),sep = " " ) %>%
          dplyr::select(V2_V3,V7, V8) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate(Control = ifelse(Genotype == "Non Infected", 1, 0),
                 Genotype = gsub("[[:punct:]]*", "", Genotype),
                 Genotype = as.character(Genotype),
                 KP = ifelse(Genotype == "KP", 1, 0),
                 Time = gsub("T", "", Time),
                 Time = as.numeric(Time))
        pd.desIb <- Ib %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>%
          as.data.frame() %>%
          dplyr::select(V8,V10, V11) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate( Control = ifelse(Genotype == "Non Infected", 1, 0),
                  Genotype = gsub("[[:punct:]]*", "", Genotype),
                  Genotype = as.character(Genotype),
                  KP = ifelse(Genotype == "KP", 1, 0),
                  Time = gsub("T", "", Time),
                  Time = as.numeric(Time))
        pd.desb<- bind_rows(pd.desNIb,pd.desIb)
        pd.desb$Replicate<- rep(1:10, 1, each = 3)
        DT::datatable(head(pd.desb), caption = "Design matrix")
        rownames(pd.desb) <- pd$geo_accession[c(4:18,34:48)]
        require(maSigPro)
        designb <- make.design.matrix(pd.desb[, c(2,3,4,5)], degree = 2)
        designb$edesign %>% as.data.frame() %>% DT::datatable()
        message("Defining the regression model...")
        designb$groups.vector
        message("sourcing functions from github...")
        require(doParallel)
        require(itertools)
        require(devtools)
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/Tfitp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/pvectorp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/position.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2sym.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2entrez.R")
        message("Fitting the regression model...")
        fitb<- p.vectorp(expb, designb,Q = 0.1, MT.adjust = "BH", min.obs = 20)
        tstepb <- T.fitp(fitb, step.method = "backward", alfa = 0.05)
        sigsb <- get.siggenes(tstepb, rsq = 0.6, vars = "groups")
        sigsb %>% glimpse()
        names(sigsb)
        names(sigsb$sig.genes)
        suma2Venn(sigsb$summary[,c(1:2)])
        clusters=see.genes(sigsb$sig.genes$KPvsControl, show.fit = T, dis=designb$dis, cluster.method = "hclust", cluster.data = 1, k.mclust = TRUE)
        df=clusters$cut %>% as.data.frame() 
        write.csv(df, "KPvsC_B_Geneclusters.csv")
        KPvsC   <- sigsb$sig.genes$KPvsControl$sig.pvalues
        message ("from probes to entrez id ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "mmusculus_gene_ensembl",
                           host = 'www.ensembl.org')
        KPvsC   <- probe2entrez(KPvsC, mart.mm)
        write.csv(KPvsC, "KPvsC_B_v0.csv")
        message ("from probes to genes symbols ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ensembl",
                         dataset = "mmusculus_gene_ensembl",
                         host="https://dec2021.archive.ensembl.org")
        mart.hs = useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                        host="https://dec2021.archive.ensembl.org")
        KPvsC   <- probe2sym(KPvsC, mart.mm,mart.hs)
        message ("writing csv files ....")
        write.csv(KPvsC, "KPvsC_B.csv")
        message ("Analysis KPvsC for Balb completed successfully ....")
      }else{
        celFiles <- list.celfiles("./gse31997data",full.names = TRUE)
        message("getting targets data from github ...")
        url <- "https://raw.githubusercontent.com/HediaTnani/Leishmania/master/targets.csv"
        targets <- getURL(url)                
        my.targets <- read.AnnotatedDataFrame(textConnection(targets),header = TRUE, row.names = 1, sep=",")
        message("Listing celfiles...")
        rawData <- read.celfiles(celFiles, phenoData = my.targets)
        colnames(rawData) <-rownames(pData(rawData))<-my.targets@data$ShortName
        eset_rma <- oligo::rma(rawData) # specify rma from package oligo or affy
        exp <- data.frame(exprs(eset_rma))
        expb <- exp[,c(4:18,34:48)]
        message("Viewing expression data ...")
        expb %>%
          as.data.frame() %>%
          head()%>%
          dplyr::select(1:4)
        message("Getting phenotypic data ...")
        gse <- getGEO("GSE31997", AnnotGPL = TRUE)
        pd <- pData(gse$GSE31997_series_matrix.txt.gz)
        NIb= as.character(pd$title[4:18])
        Ib= as.character(pd$title[34:48])
        message("Printing phenotypic data ...")
        print(NIb)
        print(Ib)
        message("Creating a design matrix ...")
        pd.desNIb <- NIb %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>%
          as.data.frame() %>%
          unite("V2_V3", c("V2","V3"),sep = " " ) %>%
          dplyr::select(V2_V3,V7, V8) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate(Control = ifelse(Genotype == "Non Infected", 1, 0),
                 Genotype = gsub("[[:punct:]]*", "", Genotype),
                 Genotype = as.character(Genotype),
                 KP = ifelse(Genotype == "KP", 1, 0),
                 Time = gsub("T", "", Time),
                 Time = as.numeric(Time))
        pd.desIb <- Ib %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>%
          as.data.frame() %>%
          dplyr::select(V8,V10, V11) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate( Control = ifelse(Genotype == "Non Infected", 1, 0),
                  Genotype = gsub("[[:punct:]]*", "", Genotype),
                  Genotype = as.character(Genotype),
                  KP = ifelse(Genotype == "KP", 1, 0),
                  Time = gsub("T", "", Time),
                  Time = as.numeric(Time))
        pd.desb<- bind_rows(pd.desNIb,pd.desIb)
        pd.desb$Replicate<- rep(1:10, 1, each = 3)
        DT::datatable(head(pd.desb), caption = "Design matrix")
        rownames(pd.desb) <- pd$geo_accession[c(4:18,34:48)]
        require(maSigPro)
        designb <- make.design.matrix(pd.desb[, c(2,3,4,5)], degree = 2)
        designb$edesign %>% as.data.frame() %>% DT::datatable()
        message("Defining the regression model...")
        designb$groups.vector
        message("sourcing functions from github...")
        require(doParallel)
        require(itertools)
        require(devtools)
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/Tfitp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/pvectorp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/position.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2sym.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2entrez.R")
        message("Fitting the regression model...")
        fitb<- p.vectorp(expb, designb,Q = 0.1, MT.adjust = "BH", min.obs = 20)
        tstepb <- T.fitp(fitb, step.method = "backward", alfa = 0.05)
        sigsb <- get.siggenes(tstepb, rsq = 0.6, vars = "groups")
        sigsb %>% glimpse()
        names(sigsb)
        names(sigsb$sig.genes)
        suma2Venn(sigsb$summary[,c(1:2)])
        clusters=see.genes(sigsb$sig.genes$KPvsControl, show.fit = T, dis=designb$dis, cluster.method = "hclust", cluster.data = 1, k.mclust = TRUE)
        df=clusters$cut %>% as.data.frame() 
        write.csv(df, "KPvsC_B_Geneclusters.csv")
        KPvsC   <- sigsb$sig.genes$KPvsControl$sig.pvalues
        message ("from probes to entrez id ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "mmusculus_gene_ensembl",
                           host = 'www.ensembl.org')
        KPvsC   <- probe2entrez(KPvsC, mart.mm)
        write.csv(KPvsC, "KPvsC_B_v0.csv")
        message ("from probes to genes symbols ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ensembl",
                         dataset = "mmusculus_gene_ensembl",
                         host="https://dec2021.archive.ensembl.org")
        mart.hs = useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                        host="https://dec2021.archive.ensembl.org")
        KPvsC   <- probe2sym(KPvsC, mart.mm, mart.hs)
        message ("writing csv files ....")
        write.csv(KPvsC, "KPvsC_B.csv")
        message ("Analysis KPvsC for Balb completed successfully ....")
      }#else
    }else if (comp== "PvsKP"){
      if (!dir.exists(paste0(getwd(),"/gse31997data"))){
        message("downloading data ...")
        getGEOSuppFiles("GSE31997")
        untar("GSE31997/GSE31997_RAW.tar", exdir="./gse31997data")
        cels <- list.files("gse31997data/", pattern = "[gz]")
        message("unzipping files ...")
        sapply(paste("gse31997data", cels, sep="/"), gunzip)
        celFiles <- list.celfiles("./gse31997data",full.names = TRUE)
        message("getting targets data from github ...")
        url <- "https://raw.githubusercontent.com/HediaTnani/Leishmania/master/targets.csv"
        targets <- getURL(url)                
        my.targets <- read.AnnotatedDataFrame(textConnection(targets),header = TRUE, row.names = 1, sep=",")
        message("Listing celfiles...")
        rawData <- read.celfiles(celFiles, phenoData = my.targets)
        colnames(rawData) <-rownames(pData(rawData))<-my.targets@data$ShortName
        eset_rma <- oligo::rma(rawData) # specify rma from package oligo or affy
        exp <- data.frame(exprs(eset_rma))
        expb <- exp[,19:48]
        message("Viewing expression data ...")
        expb %>%
          as.data.frame() %>%
          head()%>%
          dplyr::select(1:4)
        message("Getting phenotypic data ...")
        gse <- getGEO("GSE31997", AnnotGPL = TRUE)
        pd <- pData(gse$GSE31997_series_matrix.txt.gz)
        Ib= as.character(pd$title[19:33])
        NIb= as.character(pd$title[34:48])
        message("Printing phenotypic data ...")
        print(NIb)
        print(Ib)
        message("Creating a design matrix ...")
        pd.desNIb <- NIb %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>%
          as.data.frame() %>%
          dplyr::select(V8,V10, V11) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate(Genotype = gsub("[[:punct:]]*", "", Genotype),
                 Genotype = as.character(Genotype),
                 KP = ifelse(Genotype == "KP", 1, 0),
                 P = ifelse(Genotype == "P", 1, 0),
                 Time = gsub("T", "", Time),
                 Time = as.numeric(Time))
        pd.desIb <- Ib %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>%
          as.data.frame() %>%
          dplyr::select(V8,V10, V11) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate( Genotype = gsub("[[:punct:]]*", "", Genotype),
                  Genotype = as.character(Genotype),
                  KP = ifelse(Genotype == "KP", 1, 0),
                  P = ifelse(Genotype == "P", 1, 0),
                  Time = gsub("T", "", Time),
                  Time = as.numeric(Time))
        pd.desb<- bind_rows(pd.desNIb,pd.desIb)
        pd.desb$Replicate<- rep(1:10, 1, each = 3)
        DT::datatable(head(pd.desb), caption = "Design matrix")
        rownames(pd.desb) <- pd$geo_accession[19:48]
        require(maSigPro)
        designb <- make.design.matrix(pd.desb[, c(2,3,4,5)], degree = 2)
        designb$edesign %>% as.data.frame() %>% DT::datatable()
        message("Defining the regression model...")
        designb$groups.vector
        message("sourcing functions from github...")
        require(doParallel)
        require(itertools)
        require(devtools)
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/Tfitp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/pvectorp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/position.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2sym.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2entrez.R")
        message("Fitting the regression model...")
        fitb<- p.vectorp(expb, designb,Q = 0.1, MT.adjust = "BH", min.obs = 20)
        tstepb <- T.fitp(fitb, step.method = "backward", alfa = 0.05)
        sigsb <- get.siggenes(tstepb, rsq = 0.6, vars = "groups")
        sigsb %>% glimpse()
        names(sigsb)
        names(sigsb$sig.genes)
        suma2Venn(sigsb$summary[,c(1:2)])
        clusters=see.genes(sigsb$sig.genes$PvsKP, show.fit = T, dis=designb$dis, cluster.method = "hclust", cluster.data = 1, k.mclust = TRUE)
        df=clusters$cut %>% as.data.frame() 
        write.csv(df, "PvsKP_B_Geneclusters.csv")
        PvsKP <- sigsb$sig.genes$PvsKP$sig.pvalues
        message ("from probes to entrez id ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "mmusculus_gene_ensembl",
                           host = 'www.ensembl.org')
        PvsKP  <- probe2entrez(PvsKP, mart.mm)
        write.csv(PvsKP, "PvsKP_B_v0.csv")
        message ("from probes to genes symbols ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ensembl",
                         dataset = "mmusculus_gene_ensembl",
                         host="https://dec2021.archive.ensembl.org")
        mart.hs = useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                        host="https://dec2021.archive.ensembl.org")
        PvsKP   <- probe2sym(PvsKP, mart.mm, mart.hs)
        message ("writing csv files ....")
        write.csv(PvsKP, "PvsKP_B.csv")
        message ("Analysis PvsKP for Balb completed successfully ....")  
      }else{
        celFiles <- list.celfiles("./gse31997data",full.names = TRUE)
        message("getting targets data from github ...")
        url <- "https://raw.githubusercontent.com/HediaTnani/Leishmania/master/targets.csv"
        targets <- getURL(url)                
        my.targets <- read.AnnotatedDataFrame(textConnection(targets),header = TRUE, row.names = 1, sep=",")
        message("Listing celfiles...")
        rawData <- read.celfiles(celFiles, phenoData = my.targets)
        colnames(rawData) <-rownames(pData(rawData))<-my.targets@data$ShortName
        eset_rma <- oligo::rma(rawData) # specify rma from package oligo or affy
        exp <- data.frame(exprs(eset_rma))
        expb <- exp[,19:48]
        message("Viewing expression data ...")
        expb %>%
          as.data.frame() %>%
          head()%>%
          dplyr::select(1:4)
        message("Getting phenotypic data ...")
        gse <- getGEO("GSE31997", AnnotGPL = TRUE)
        pd <- pData(gse$GSE31997_series_matrix.txt.gz)
        Ib= as.character(pd$title[19:33])
        NIb= as.character(pd$title[34:48])
        message("Printing phenotypic data ...")
        print(NIb)
        print(Ib)
        message("Creating a design matrix ...")
        pd.desNIb <- NIb %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>%
          as.data.frame() %>%
          dplyr::select(V8,V10, V11) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate(Genotype = gsub("[[:punct:]]*", "", Genotype),
                 Genotype = as.character(Genotype),
                 KP = ifelse(Genotype == "KP", 1, 0),
                 P = ifelse(Genotype == "P", 1, 0),
                 Time = gsub("T", "", Time),
                 Time = as.numeric(Time))
        pd.desIb <- Ib %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>%
          as.data.frame() %>%
          dplyr::select(V8,V10, V11) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate( Genotype = gsub("[[:punct:]]*", "", Genotype),
                  Genotype = as.character(Genotype),
                  KP = ifelse(Genotype == "KP", 1, 0),
                  P = ifelse(Genotype == "P", 1, 0),
                  Time = gsub("T", "", Time),
                  Time = as.numeric(Time))
        pd.desb<- bind_rows(pd.desNIb,pd.desIb)
        pd.desb$Replicate<- rep(1:10, 1, each = 3)
        DT::datatable(head(pd.desb), caption = "Design matrix")
        rownames(pd.desb) <- pd$geo_accession[19:48]
        require(maSigPro)
        designb <- make.design.matrix(pd.desb[, c(2,3,4,5)], degree = 2)
        designb$edesign %>% as.data.frame() %>% DT::datatable()
        message("Defining the regression model...")
        designb$groups.vector
        message("sourcing functions from github...")
        require(doParallel)
        require(itertools)
        require(devtools)
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/Tfitp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/pvectorp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/position.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2sym.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2entrez.R")
        message("Fitting the regression model...")
        fitb<- p.vectorp(expb, designb,Q = 0.1, MT.adjust = "BH", min.obs = 20)
        tstepb <- T.fitp(fitb, step.method = "backward", alfa = 0.05)
        sigsb <- get.siggenes(tstepb, rsq = 0.6, vars = "groups")
        sigsb %>% glimpse()
        names(sigsb)
        names(sigsb$sig.genes)
        suma2Venn(sigsb$summary[,c(1:2)])
        clusters=see.genes(sigsb$sig.genes$PvsKP, show.fit = T, dis=designb$dis, cluster.method = "hclust", cluster.data = 1, k.mclust = TRUE)
        df=clusters$cut %>% as.data.frame() 
        write.csv(df, "PvsKP_B_Geneclusters.csv")
        PvsKP   <- sigsb$sig.genes$PvsKP$sig.pvalues
        message ("from probes to entrez id ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "mmusculus_gene_ensembl",
                           host = 'www.ensembl.org')
        PvsKP   <- probe2entrez(PvsKP, mart.mm)
        write.csv(PvsKP, "PvsKP_B_v0.csv")
        message ("from probes to genes symbols ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ensembl",
                         dataset = "mmusculus_gene_ensembl",
                         host="https://dec2021.archive.ensembl.org")
        mart.hs = useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                        host="https://dec2021.archive.ensembl.org")
        PvsKP   <- probe2sym(PvsKP, mart.mm, mart.hs)
        message ("writing csv files ....")
        write.csv(PvsKP, "PvsKP_B.csv")
        message ("Analysis PvsKP for Balb completed successfully ....")
      }#else
    }#PvsKP
  }#Balb
  if(mouse=="c57"){#c57
    message ("loading packages ....")
    require(GEOquery)
    require(oligo)
    require(mogene10sttranscriptcluster.db)
    require(RCurl)
    require(foreign)
    require(tidyverse)
    setwd(getwd())
    if (comp== "IvsC"){
      if (!dir.exists(paste0(getwd(),"/gse31997data"))){
        message("downloading data ...")
        getGEOSuppFiles("GSE31997")
        untar("GSE31997/GSE31997_RAW.tar", exdir="./gse31997data")
        cels <- list.files("gse31997data/", pattern = "[gz]")
        message("unzipping files ...")
        sapply(paste("gse31997data", cels, sep="/"), gunzip)
        celFiles <- list.celfiles("./gse31997data",full.names = TRUE)
        message("getting targets data from github ...")
        url <- "https://raw.githubusercontent.com/HediaTnani/Leishmania/master/targets.csv"
        targets <- getURL(url)                
        my.targets <- read.AnnotatedDataFrame(textConnection(targets),header = TRUE, row.names = 1, sep=",")
        message("Listing celfiles...")
        rawData <- read.celfiles(celFiles, phenoData = my.targets)
        colnames(rawData) <-rownames(pData(rawData))<-my.targets@data$ShortName
        eset_rma <- oligo::rma(rawData) # specify rma from package oligo or affy
        exp <- data.frame(exprs(eset_rma))
        expc <- exp[,52:96]
        message("Viewing expression data ...")
        expc %>%
          as.data.frame() %>%
          head()%>%
          dplyr::select(1:4) 
        message("Getting phenotypic data ...")
        gse <- getGEO("GSE31997", AnnotGPL = TRUE)
        pd <- pData(gse$GSE31997_series_matrix.txt.gz)
        NIc= as.character(pd$title[52:66]) 
        Ic= as.character(pd$title[67:96]) 
        message("Printing phenotypic data ...")
        print(NIc)
        print(Ic)
        message("Creating a design matrix ...")
        pd.desNIc <- NIc %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
          as.data.frame() %>% 
          unite("V2_V3", c("V2","V3"),sep = " " ) %>%
          dplyr::select(V2_V3,V7, V8) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate(Control = ifelse(Genotype == "Non Infected", 1, 0),
                 Genotype = gsub("[[:punct:]]*", "", Genotype),
                 Genotype = as.character(Genotype),
                 Infected = ifelse(Genotype == "P"|Genotype == "KP", 1, 0),
                 Time = gsub("T", "", Time),
                 Time = as.numeric(Time))
        pd.desIc <- Ic %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
          as.data.frame() %>% 
          dplyr::select(V8,V10, V11) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate( Control = ifelse(Genotype == "Non Infected", 1, 0),
                  Genotype = gsub("[[:punct:]]*", "", Genotype),
                  Genotype = as.character(Genotype),
                  Infected = ifelse(Genotype == "P"|Genotype == "KP", 1, 0),
                  Time = gsub("T", "", Time),
                  Time = as.numeric(Time))
        pd.desc<- bind_rows(pd.desNIc,pd.desIc)
        pd.desc$Replicate<- rep(1:15, 1, each = 3)
        DT::datatable(head(pd.desc), caption = "Design matrix")
        rownames(pd.desc) <- pd$geo_accession[52:96]
        require(maSigPro)
        designc <- make.design.matrix(pd.desc[, c(2,3,4,5)], degree = 2)
        designc$edesign %>%
          as.data.frame() %>% 
          DT::datatable()
        message("Defining the regression model...")
        designc$groups.vector
        message("sourcing functions from github...")
        require(doParallel)
        require(itertools)
        require(devtools)
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/Tfitp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/pvectorp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/position.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2sym.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2entrez.R")
        message("Fitting the regression model...")
        fitc<- p.vectorp(expc, designc,Q = 0.1, MT.adjust = "BH", min.obs = 3)
        tstepc <- T.fitp(fitc, step.method = "backward", alfa = 0.05)
        sigsc <- get.siggenes(tstepc, rsq = 0.6, vars = "groups")
        sigsc %>%
          glimpse()
        names(sigsc)
        names(sigsc$sig.genes)
        suma2Venn(sigsc$summary[,c(1:2)])
        clusters=see.genes(sigsc$sig.genes$InfectedvsControl, show.fit = T, dis=designc$dis, cluster.method = "hclust", cluster.data = 1, k.mclust = TRUE)
        df=clusters$cut %>% as.data.frame() 
        write.csv(df, "IvsC_C_Geneclusters.csv")
        IvsC   <- sigsc$sig.genes$InfectedvsControl$sig.pvalues
        message ("from probes to entrez id ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "mmusculus_gene_ensembl",
                           host = 'www.ensembl.org')
        IvsC   <- probe2entrez(IvsC, mart.mm)
        write.csv(IvsC, "IvsC_C_v0.csv")
        message ("from probes to genes symbols ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ensembl",
                         dataset = "mmusculus_gene_ensembl",
                         host="https://dec2021.archive.ensembl.org")
        mart.hs = useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                        host="https://dec2021.archive.ensembl.org")
        IvsC   <- probe2sym(IvsC, mart.mm, mart.hs)
        message ("writing csv files ....")
        write.csv(IvsC, "IvsC_C.csv")
        message ("Analysis IvsC for c57 completed successfully ....")
      }else{
        celFiles <- list.celfiles("./gse31997data",full.names = TRUE)
        message("getting targets data from github ...")
        url <- "https://raw.githubusercontent.com/HediaTnani/Leishmania/master/targets.csv"
        targets <- getURL(url)                
        my.targets <- read.AnnotatedDataFrame(textConnection(targets),header = TRUE, row.names = 1, sep=",")
        message("Listing celfiles...")
        rawData <- read.celfiles(celFiles, phenoData = my.targets)
        colnames(rawData) <-rownames(pData(rawData))<-my.targets@data$ShortName
        eset_rma <- oligo::rma(rawData) # specify rma from package oligo or affy
        exp <- data.frame(exprs(eset_rma))
        expc <- exp[,52:96]
        message("Viewing expression data ...")
        expc %>%
          as.data.frame() %>%
          head()%>%
          dplyr::select(1:4) 
        message("Getting phenotypic data ...")
        gse <- getGEO("GSE31997", AnnotGPL = TRUE)
        pd <- pData(gse$GSE31997_series_matrix.txt.gz)
        NIc= as.character(pd$title[52:66]) 
        Ic= as.character(pd$title[67:96]) 
        message("Printing phenotypic data ...")
        print(NIc)
        print(Ic)
        message("Creating a design matrix ...")
        pd.desNIc <- NIc %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
          as.data.frame() %>% 
          unite("V2_V3", c("V2","V3"),sep = " " ) %>%
          dplyr::select(V2_V3,V7, V8) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate(Control = ifelse(Genotype == "Non Infected", 1, 0),
                 Genotype = gsub("[[:punct:]]*", "", Genotype),
                 Genotype = as.character(Genotype),
                 Infected = ifelse(Genotype == "P"|Genotype == "KP", 1, 0),
                 Time = gsub("T", "", Time),
                 Time = as.numeric(Time))
        pd.desIc <- Ic %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
          as.data.frame() %>% 
          dplyr::select(V8,V10, V11) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate( Control = ifelse(Genotype == "Non Infected", 1, 0),
                  Genotype = gsub("[[:punct:]]*", "", Genotype),
                  Genotype = as.character(Genotype),
                  Infected = ifelse(Genotype == "P"|Genotype == "KP", 1, 0),
                  Time = gsub("T", "", Time),
                  Time = as.numeric(Time))
        pd.desc<- bind_rows(pd.desNIc,pd.desIc)
        pd.desc$Replicate<- rep(1:15, 1, each = 3)
        DT::datatable(head(pd.desc), caption = "Design matrix")
        rownames(pd.desc) <- pd$geo_accession[52:96]
        require(maSigPro)
        designc <- make.design.matrix(pd.desc[, c(2,3,4,5)], degree = 2)
        designc$edesign %>%
          as.data.frame() %>% 
          DT::datatable()
        message("Defining the regression model...")
        designc$groups.vector
        message("sourcing functions from github...")
        require(doParallel)
        require(itertools)
        require(devtools)
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/Tfitp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/pvectorp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/position.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2sym.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2entrez.R")
        message("Fitting the regression model...")
        fitc<- p.vectorp(expc, designc,Q = 0.1, MT.adjust = "BH", min.obs = 3)
        tstepc <- T.fitp(fitc, step.method = "backward", alfa = 0.05)
        sigsc <- get.siggenes(tstepc, rsq = 0.6, vars = "groups")
        sigsc %>%
          glimpse()
        names(sigsc)
        names(sigsc$sig.genes)
        suma2Venn(sigsc$summary[,c(1:2)])
        clusters=see.genes(sigsc$sig.genes$InfectedvsControl, show.fit = T, dis=designc$dis, cluster.method = "hclust", cluster.data = 1, k.mclust = TRUE)
        df=clusters$cut %>% as.data.frame() 
        write.csv(df, "IvsC_C_Geneclusters.csv")
        IvsC   <- sigsc$sig.genes$InfectedvsControl$sig.pvalues
        message ("from probes to entrez id ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "mmusculus_gene_ensembl",
                           host = 'www.ensembl.org')
        IvsC   <- probe2entrez(IvsC, mart.mm)
        write.csv(IvsC, "IvsC_C_v0.csv")
        message ("from probes to genes symbols ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ensembl",
                         dataset = "mmusculus_gene_ensembl",
                         host="https://dec2021.archive.ensembl.org")
        mart.hs = useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                        host="https://dec2021.archive.ensembl.org")
        IvsC   <- probe2sym(IvsC, mart.mm, mart.hs)
        message ("writing csv files ....")
        write.csv(IvsC, "IvsC_C.csv")
        message ("Analysis IvsC for c57 completed successfully ....")
      }#else
    }else if (comp== "PvsC"){
      if (!dir.exists(paste0(getwd(),"/gse31997data"))){
        message("downloading data ...")
        getGEOSuppFiles("GSE31997")
        untar("GSE31997/GSE31997_RAW.tar", exdir="./gse31997data")
        cels <- list.files("gse31997data/", pattern = "[gz]")
        message("unzipping files ...")
        sapply(paste("gse31997data", cels, sep="/"), gunzip)
        celFiles <- list.celfiles("./gse31997data",full.names = TRUE)
        message("getting targets data from github ...")
        url <- "https://raw.githubusercontent.com/HediaTnani/Leishmania/master/targets.csv"
        targets <- getURL(url)                
        my.targets <- read.AnnotatedDataFrame(textConnection(targets),header = TRUE, row.names = 1, sep=",")
        message("Listing celfiles...")
        rawData <- read.celfiles(celFiles, phenoData = my.targets)
        colnames(rawData) <-rownames(pData(rawData))<-my.targets@data$ShortName
        eset_rma <- oligo::rma(rawData) # specify rma from package oligo or affy
        exp <- data.frame(exprs(eset_rma))
        expc <- exp[,52:81]
        message("Viewing expression data ...")
        expc %>%
          as.data.frame() %>%
          head()%>%
          dplyr::select(1:4) 
        message("Getting phenotypic data ...")
        gse <- getGEO("GSE31997", AnnotGPL = TRUE)
        pd <- pData(gse$GSE31997_series_matrix.txt.gz)
        NIc= as.character(pd$title[52:66]) 
        Ic= as.character(pd$title[67:81]) 
        message("Printing phenotypic data ...")
        print(NIc)
        print(Ic)
        message("Creating a design matrix ...")
        pd.desNIc <- NIc %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
          as.data.frame() %>% 
          unite("V2_V3", c("V2","V3"),sep = " " ) %>%
          dplyr::select(V2_V3,V7, V8) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate(Control = ifelse(Genotype == "Non Infected", 1, 0),
                 Genotype = gsub("[[:punct:]]*", "", Genotype),
                 Genotype = as.character(Genotype),
                 P = ifelse(Genotype == "P", 1, 0),
                 Time = gsub("T", "", Time),
                 Time = as.numeric(Time))
        pd.desIc <- Ic %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
          as.data.frame() %>% 
          dplyr::select(V8,V10, V11) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate( Control = ifelse(Genotype == "Non Infected", 1, 0),
                  Genotype = gsub("[[:punct:]]*", "", Genotype),
                  Genotype = as.character(Genotype),
                  P = ifelse(Genotype == "P", 1, 0),
                  Time = gsub("T", "", Time),
                  Time = as.numeric(Time))
        pd.desc<- bind_rows(pd.desNIc,pd.desIc)
        pd.desc$Replicate<- rep(1:10, 1, each = 3)
        DT::datatable(head(pd.desc), caption = "Design matrix")
        rownames(pd.desc) <- pd$geo_accession[52:81]
        require(maSigPro)
        designc <- make.design.matrix(pd.desc[, c(2,3,4,5)], degree = 2)
        designc$edesign %>%
          as.data.frame() %>% 
          DT::datatable()
        message("Defining the regression model...")
        designc$groups.vector
        message("sourcing functions from github...")
        require(doParallel)
        require(itertools)
        require(devtools)
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/Tfitp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/pvectorp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/position.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2sym.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2entrez.R")
        message("Fitting the regression model...")
        fitc<- p.vectorp(expc, designc,Q = 0.1, MT.adjust = "BH", min.obs = 3)
        tstepc <- T.fitp(fitc, step.method = "backward", alfa = 0.05)
        sigsc <- get.siggenes(tstepc, rsq = 0.6, vars = "groups")
        sigsc %>%
          glimpse()
        names(sigsc)
        names(sigsc$sig.genes)
        suma2Venn(sigsc$summary[,c(1:2)])
        clusters=see.genes(sigsc$sig.genes$PvsControl, show.fit = T, dis=designc$dis, cluster.method = "hclust", cluster.data = 1, k.mclust = TRUE)
        df=clusters$cut %>% as.data.frame() 
        write.csv(df, "PvsC_C_Geneclusters.csv")
        PvsC   <- sigsc$sig.genes$PvsControl$sig.pvalues
        message ("from probes to entrez id ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "mmusculus_gene_ensembl",
                           host = 'www.ensembl.org')
        PvsC   <- probe2entrez(PvsC, mart.mm)
        write.csv(PvsC, "PvsC_C_v0.csv")
        message ("from probes to genes symbols ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ensembl",
                         dataset = "mmusculus_gene_ensembl",
                         host="https://dec2021.archive.ensembl.org")
        mart.hs = useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                        host="https://dec2021.archive.ensembl.org")
        PvsC   <- probe2sym(PvsC, mart.mm, mart.hs)
        message ("writing csv files ....")
        write.csv(PvsC, "PvsC_C.csv")
        message ("Analysis PvsC for c57 completed successfully ....")
      }else{
        celFiles <- list.celfiles("./gse31997data",full.names = TRUE)
        message("getting targets data from github ...")
        url <- "https://raw.githubusercontent.com/HediaTnani/Leishmania/master/targets.csv"
        targets <- getURL(url)                
        my.targets <- read.AnnotatedDataFrame(textConnection(targets),header = TRUE, row.names = 1, sep=",")
        message("Listing celfiles...")
        rawData <- read.celfiles(celFiles, phenoData = my.targets)
        colnames(rawData) <-rownames(pData(rawData))<-my.targets@data$ShortName
        eset_rma <- oligo::rma(rawData) # specify rma from package oligo or affy
        exp <- data.frame(exprs(eset_rma))
        expc <- exp[,52:81]
        message("Viewing expression data ...")
        expc %>%
          as.data.frame() %>%
          head()%>%
          dplyr::select(1:4) 
        message("Getting phenotypic data ...")
        gse <- getGEO("GSE31997", AnnotGPL = TRUE)
        pd <- pData(gse$GSE31997_series_matrix.txt.gz)
        NIc= as.character(pd$title[52:66]) 
        Ic= as.character(pd$title[67:81]) 
        message("Printing phenotypic data ...")
        print(NIc)
        print(Ic)
        message("Creating a design matrix ...")
        pd.desNIc <- NIc %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
          as.data.frame() %>% 
          unite("V2_V3", c("V2","V3"),sep = " " ) %>%
          dplyr::select(V2_V3,V7, V8) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate(Control = ifelse(Genotype == "Non Infected", 1, 0),
                 Genotype = gsub("[[:punct:]]*", "", Genotype),
                 Genotype = as.character(Genotype),
                 P = ifelse(Genotype == "P", 1, 0),
                 Time = gsub("T", "", Time),
                 Time = as.numeric(Time))
        pd.desIc <- Ic %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
          as.data.frame() %>% 
          dplyr::select(V8,V10, V11) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate( Control = ifelse(Genotype == "Non Infected", 1, 0),
                  Genotype = gsub("[[:punct:]]*", "", Genotype),
                  Genotype = as.character(Genotype),
                  P = ifelse(Genotype == "P", 1, 0),
                  Time = gsub("T", "", Time),
                  Time = as.numeric(Time))
        pd.desc<- bind_rows(pd.desNIc,pd.desIc)
        pd.desc$Replicate<- rep(1:10, 1, each = 3)
        DT::datatable(head(pd.desc), caption = "Design matrix")
        rownames(pd.desc) <- pd$geo_accession[52:81]
        require(maSigPro)
        designc <- make.design.matrix(pd.desc[, c(2,3,4,5)], degree = 2)
        designc$edesign %>%
          as.data.frame() %>% 
          DT::datatable()
        message("Defining the regression model...")
        designc$groups.vector
        message("sourcing functions from github...")
        require(doParallel)
        require(itertools)
        require(devtools)
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/Tfitp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/pvectorp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/position.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2sym.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2entrez.R")
        message("Fitting the regression model...")
        fitc<- p.vectorp(expc, designc,Q = 0.1, MT.adjust = "BH", min.obs = 3)
        tstepc <- T.fitp(fitc, step.method = "backward", alfa = 0.05)
        sigsc <- get.siggenes(tstepc, rsq = 0.6, vars = "groups")
        sigsc %>%
          glimpse()
        names(sigsc)
        names(sigsc$sig.genes)
        suma2Venn(sigsc$summary[,c(1:2)])
        clusters=see.genes(sigsc$sig.genes$PvsControl, show.fit = T, dis=designc$dis, cluster.method = "hclust", cluster.data = 1, k.mclust = TRUE)
        df=clusters$cut %>% as.data.frame() 
        write.csv(df, "PvsC_C_Geneclusters.csv")
        PvsC   <- sigsc$sig.genes$PvsControl$sig.pvalues
        message ("from probes to entrez id ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "mmusculus_gene_ensembl",
                           host = 'www.ensembl.org')
        PvsC   <- probe2entrez(PvsC, mart.mm)
        write.csv(PvsC, "PvsC_C_v0.csv")
        message ("from probes to genes symbols ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ensembl",
                         dataset = "mmusculus_gene_ensembl",
                         host="https://dec2021.archive.ensembl.org")
        mart.hs = useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                        host="https://dec2021.archive.ensembl.org")
        PvsC   <- probe2sym(PvsC, mart.mm, mart.hs)
        message ("writing csv files ....")
        write.csv(PvsC, "PvsC_C.csv")
        message ("Analysis PvsC for c57 completed successfully ....")
      }#else
    }else if (comp== "KPvsC"){
      if (!dir.exists(paste0(getwd(),"/gse31997data"))){
        message("downloading data ...")
        getGEOSuppFiles("GSE31997")
        untar("GSE31997/GSE31997_RAW.tar", exdir="./gse31997data")
        cels <- list.files("gse31997data/", pattern = "[gz]")
        message("unzipping files ...")
        sapply(paste("gse31997data", cels, sep="/"), gunzip)
        celFiles <- list.celfiles("./gse31997data",full.names = TRUE)
        message("getting targets data from github ...")
        url <- "https://raw.githubusercontent.com/HediaTnani/Leishmania/master/targets.csv"
        targets <- getURL(url)                
        my.targets <- read.AnnotatedDataFrame(textConnection(targets),header = TRUE, row.names = 1, sep=",")
        message("Listing celfiles...")
        rawData <- read.celfiles(celFiles, phenoData = my.targets)
        colnames(rawData) <-rownames(pData(rawData))<-my.targets@data$ShortName
        eset_rma <- oligo::rma(rawData) # specify rma from package oligo or affy
        exp <- data.frame(exprs(eset_rma))
        expc <- exp[,c(52:66,82:96)]
        message("Viewing expression data ...")
        expc %>%
          as.data.frame() %>%
          head()%>%
          dplyr::select(1:4) 
        message("Getting phenotypic data ...")
        gse <- getGEO("GSE31997", AnnotGPL = TRUE)
        pd <- pData(gse$GSE31997_series_matrix.txt.gz)
        NIc= as.character(pd$title[52:66]) 
        Ic= as.character(pd$title[82:96]) 
        message("Printing phenotypic data ...")
        print(NIc)
        print(Ic)
        message("Creating a design matrix ...")
        pd.desNIc <- NIc %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
          as.data.frame() %>% 
          unite("V2_V3", c("V2","V3"),sep = " " ) %>%
          dplyr::select(V2_V3,V7, V8) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate(Control = ifelse(Genotype == "Non Infected", 1, 0),
                 Genotype = gsub("[[:punct:]]*", "", Genotype),
                 Genotype = as.character(Genotype),
                 KP = ifelse(Genotype == "KP", 1, 0),
                 Time = gsub("T", "", Time),
                 Time = as.numeric(Time))
        pd.desIc <- Ic %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
          as.data.frame() %>% 
          dplyr::select(V8,V10, V11) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate( Control = ifelse(Genotype == "Non Infected", 1, 0),
                  Genotype = gsub("[[:punct:]]*", "", Genotype),
                  Genotype = as.character(Genotype),
                  KP = ifelse(Genotype == "KP", 1, 0),
                  Time = gsub("T", "", Time),
                  Time = as.numeric(Time))
        pd.desc<- bind_rows(pd.desNIc,pd.desIc)
        pd.desc$Replicate<- rep(1:10, 1, each = 3)
        DT::datatable(head(pd.desc), caption = "Design matrix")
        rownames(pd.desc) <- pd$geo_accession[c(52:66,82:96)]
        require(maSigPro)
        designc <- make.design.matrix(pd.desc[, c(2,3,4,5)], degree = 2)
        designc$edesign %>%
          as.data.frame() %>% 
          DT::datatable()
        message("Defining the regression model...")
        designc$groups.vector
        message("sourcing functions from github...")
        require(doParallel)
        require(itertools)
        require(devtools)
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/Tfitp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/pvectorp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/position.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2sym.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2entrez.R")
        message("Fitting the regression model...")
        fitc<- p.vectorp(expc, designc,Q = 0.1, MT.adjust = "BH", min.obs = 3)
        tstepc <- T.fitp(fitc, step.method = "backward", alfa = 0.05)
        sigsc <- get.siggenes(tstepc, rsq = 0.6, vars = "groups")
        sigsc %>%
          glimpse()
        names(sigsc)
        names(sigsc$sig.genes)
        suma2Venn(sigsc$summary[,c(1:2)])
        clusters=see.genes(sigsc$sig.genes$KPvsControl, show.fit = T, dis=designc$dis, cluster.method = "hclust", cluster.data = 1, k.mclust = TRUE)
        df=clusters$cut %>% as.data.frame() 
        write.csv(df, "KPvsC_C_Geneclusters.csv")
        KPvsC <- sigsc$sig.genes$KPvsControl$sig.pvalues
        message ("from probes to entrez id ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "mmusculus_gene_ensembl",
                           host = 'www.ensembl.org')
        KPvsC   <- probe2entrez(KPvsC, mart.mm)
        write.csv(KPvsC, "KPvsC_C_v0.csv")
        message ("from probes to genes symbols ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ensembl",
                         dataset = "mmusculus_gene_ensembl",
                         host="https://dec2021.archive.ensembl.org")
        mart.hs = useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                        host="https://dec2021.archive.ensembl.org")
        KPvsC <- probe2sym(KPvsC, mart.mm, mart.hs)
        message ("writing csv files ....")
        write.csv(KPvsC, "KPvsC_C.csv")
        message ("Analysis KPvsC for c57 completed successfully ....")
      }else{
        celFiles <- list.celfiles("./gse31997data",full.names = TRUE)
        message("getting targets data from github ...")
        url <- "https://raw.githubusercontent.com/HediaTnani/Leishmania/master/targets.csv"
        targets <- getURL(url)                
        my.targets <- read.AnnotatedDataFrame(textConnection(targets),header = TRUE, row.names = 1, sep=",")
        message("Listing celfiles...")
        rawData <- read.celfiles(celFiles, phenoData = my.targets)
        colnames(rawData) <-rownames(pData(rawData))<-my.targets@data$ShortName
        eset_rma <- oligo::rma(rawData) # specify rma from package oligo or affy
        exp <- data.frame(exprs(eset_rma))
        expc <- exp[,c(52:66,82:96)]
        message("Viewing expression data ...")
        expc %>%
          as.data.frame() %>%
          head()%>%
          dplyr::select(1:4) 
        message("Getting phenotypic data ...")
        gse <- getGEO("GSE31997", AnnotGPL = TRUE)
        pd <- pData(gse$GSE31997_series_matrix.txt.gz)
        NIc= as.character(pd$title[52:66]) 
        Ic= as.character(pd$title[82:96]) 
        message("Printing phenotypic data ...")
        print(NIc)
        print(Ic)
        message("Creating a design matrix ...")
        pd.desNIc <- NIc %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
          as.data.frame() %>% 
          unite("V2_V3", c("V2","V3"),sep = " " ) %>%
          dplyr::select(V2_V3,V7, V8) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate(Control = ifelse(Genotype == "Non Infected", 1, 0),
                 Genotype = gsub("[[:punct:]]*", "", Genotype),
                 Genotype = as.character(Genotype),
                 KP = ifelse(Genotype == "KP", 1, 0),
                 Time = gsub("T", "", Time),
                 Time = as.numeric(Time))
        pd.desIc <- Ic %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
          as.data.frame() %>% 
          dplyr::select(V8,V10, V11) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate( Control = ifelse(Genotype == "Non Infected", 1, 0),
                  Genotype = gsub("[[:punct:]]*", "", Genotype),
                  Genotype = as.character(Genotype),
                  KP = ifelse(Genotype == "KP", 1, 0),
                  Time = gsub("T", "", Time),
                  Time = as.numeric(Time))
        pd.desc<- bind_rows(pd.desNIc,pd.desIc)
        pd.desc$Replicate<- rep(1:10, 1, each = 3)
        DT::datatable(head(pd.desc), caption = "Design matrix")
        rownames(pd.desc) <- pd$geo_accession[c(52:66,82:96)]
        require(maSigPro)
        designc <- make.design.matrix(pd.desc[, c(2,3,4,5)], degree = 2)
        designc$edesign %>%
          as.data.frame() %>% 
          DT::datatable()
        message("Defining the regression model...")
        designc$groups.vector
        message("sourcing functions from github...")
        require(doParallel)
        require(itertools)
        require(devtools)
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/Tfitp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/pvectorp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/position.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2sym.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2entrez.R")
        message("Fitting the regression model...")
        fitc<- p.vectorp(expc, designc,Q = 0.1, MT.adjust = "BH", min.obs = 3)
        tstepc <- T.fitp(fitc, step.method = "backward", alfa = 0.05)
        sigsc <- get.siggenes(tstepc, rsq = 0.6, vars = "groups")
        sigsc %>%
          glimpse()
        names(sigsc)
        names(sigsc$sig.genes)
        suma2Venn(sigsc$summary[,c(1:2)])
        clusters=see.genes(sigsc$sig.genes$KPvsControl, show.fit = T, dis=designc$dis, cluster.method = "hclust", cluster.data = 1, k.mclust = TRUE)
        df=clusters$cut %>% as.data.frame() 
        write.csv(df, "KPvsC_C_Geneclusters.csv")
        KPvsC <- sigsc$sig.genes$KPvsControl$sig.pvalues
        message ("from probes to entrez id ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "mmusculus_gene_ensembl",
                           host = 'www.ensembl.org')
        KPvsC   <- probe2entrez(KPvsC, mart.mm)
        write.csv(KPvsC, "KPvsC_C_v0.csv")
        message ("from probes to genes symbols ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ensembl",
                         dataset = "mmusculus_gene_ensembl",
                         host="https://dec2021.archive.ensembl.org")
        mart.hs = useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                        host="https://dec2021.archive.ensembl.org")
        KPvsC <- probe2sym(KPvsC, mart.mm, mart.hs)
        message ("writing csv files ....")
        write.csv(KPvsC, "KPvsC_C.csv")
        message ("Analysis KPvsC for C57 completed successfully ....")
      }#else
    }else if (comp== "PvsKP"){
      if (!dir.exists(paste0(getwd(),"/gse31997data"))){
        message("downloading data ...")
        getGEOSuppFiles("GSE31997")
        untar("GSE31997/GSE31997_RAW.tar", exdir="./gse31997data")
        cels <- list.files("gse31997data/", pattern = "[gz]")
        message("unzipping files ...")
        sapply(paste("gse31997data", cels, sep="/"), gunzip)
        celFiles <- list.celfiles("./gse31997data",full.names = TRUE)
        message("getting targets data from github ...")
        url <- "https://raw.githubusercontent.com/HediaTnani/Leishmania/master/targets.csv"
        targets <- getURL(url)                
        my.targets <- read.AnnotatedDataFrame(textConnection(targets),header = TRUE, row.names = 1, sep=",")
        message("Listing celfiles...")
        rawData <- read.celfiles(celFiles, phenoData = my.targets)
        colnames(rawData) <-rownames(pData(rawData))<-my.targets@data$ShortName
        eset_rma <- oligo::rma(rawData) # specify rma from package oligo or affy
        exp <- data.frame(exprs(eset_rma))
        expc <- exp[,67:96]
        message("Viewing expression data ...")
        expc %>%
          as.data.frame() %>%
          head()%>%
          dplyr::select(1:4) 
        message("Getting phenotypic data ...")
        gse <- getGEO("GSE31997", AnnotGPL = TRUE)
        pd <- pData(gse$GSE31997_series_matrix.txt.gz)
        NIc= as.character(pd$title[67:81]) 
        Ic= as.character(pd$title[82:96]) 
        message("Printing phenotypic data ...")
        print(NIc)
        print(Ic)
        message("Creating a design matrix ...")
        pd.desNIc <- NIc %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
          as.data.frame() %>% 
          dplyr::select(V8,V10, V11) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate( Genotype = gsub("[[:punct:]]*", "", Genotype),
                  Genotype = as.character(Genotype),
                  KP = ifelse(Genotype == "KP", 1, 0),
                  P = ifelse(Genotype == "P", 1, 0),
                  Time = gsub("T", "", Time),
                  Time = as.numeric(Time))
        pd.desIc <- Ic %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
          as.data.frame() %>% 
          dplyr::select(V8,V10, V11) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate( Genotype = gsub("[[:punct:]]*", "", Genotype),
                  Genotype = as.character(Genotype),
                  KP = ifelse(Genotype == "KP", 1, 0),
                  P = ifelse(Genotype == "P", 1, 0),
                  Time = gsub("T", "", Time),
                  Time = as.numeric(Time))
        pd.desc<- bind_rows(pd.desNIc,pd.desIc)
        pd.desc$Replicate<- rep(1:10, 1, each = 3)
        DT::datatable(head(pd.desc), caption = "Design matrix")
        rownames(pd.desc) <- pd$geo_accession[67:96]
        require(maSigPro)
        designc <- make.design.matrix(pd.desc[, c(2,3,4,5)], degree = 2)
        designc$edesign %>%
          as.data.frame() %>% 
          DT::datatable()
        message("Defining the regression model...")
        designc$groups.vector
        message("sourcing functions from github...")
        require(doParallel)
        require(itertools)
        require(devtools)
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/Tfitp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/pvectorp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/position.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2sym.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2entrez.R")
        message("Fitting the regression model...")
        fitc<- p.vectorp(expc, designc,Q = 0.1, MT.adjust = "BH", min.obs = 3)
        tstepc <- T.fitp(fitc, step.method = "backward", alfa = 0.05)
        sigsc <- get.siggenes(tstepc, rsq = 0.6, vars = "groups")
        sigsc %>%
          glimpse()
        names(sigsc)
        names(sigsc$sig.genes)
        suma2Venn(sigsc$summary[,c(1:2)])
        clusters=see.genes(sigsc$sig.genes$PvsKP, show.fit = T, dis=designc$dis, cluster.method = "hclust", cluster.data = 1, k.mclust = TRUE)
        df=clusters$cut %>% as.data.frame() 
        write.csv(df, "PvsKP_C_Geneclusters.csv")
        PvsKP <- sigsc$sig.genes$PvsKP$sig.pvalues
        message ("from probes to entrez id ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "mmusculus_gene_ensembl",
                           host = 'www.ensembl.org')
        PvsKP   <- probe2entrez(PvsKP, mart.mm)
        write.csv(PvsKP, "PvsKP_C_v0.csv")
        message ("from probes to genes symbols ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ensembl",
                         dataset = "mmusculus_gene_ensembl",
                         host="https://dec2021.archive.ensembl.org")
        mart.hs = useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                        host="https://dec2021.archive.ensembl.org")
        PvsKP <- probe2sym(PvsKP, mart.mm, mart.hs)
        message ("writing csv files ....")
        write.csv(PvsKP, "PvsKP_C.csv")
        message ("Analysis PvsKP completed successfully ....")
      }else{
        celFiles <- list.celfiles("./gse31997data",full.names = TRUE)
        message("getting targets data from github ...")
        url <- "https://raw.githubusercontent.com/HediaTnani/Leishmania/master/targets.csv"
        targets <- getURL(url)                
        my.targets <- read.AnnotatedDataFrame(textConnection(targets),header = TRUE, row.names = 1, sep=",")
        message("Listing celfiles...")
        rawData <- read.celfiles(celFiles, phenoData = my.targets)
        colnames(rawData) <-rownames(pData(rawData))<-my.targets@data$ShortName
        eset_rma <- oligo::rma(rawData) # specify rma from package oligo or affy
        exp <- data.frame(exprs(eset_rma))
        expc <- exp[,67:96]
        message("Viewing expression data ...")
        expc %>%
          as.data.frame() %>%
          head()%>%
          dplyr::select(1:4) 
        message("Getting phenotypic data ...")
        gse <- getGEO("GSE31997", AnnotGPL = TRUE)
        pd <- pData(gse$GSE31997_series_matrix.txt.gz)
        NIc= as.character(pd$title[67:81]) 
        Ic= as.character(pd$title[82:96]) 
        message("Printing phenotypic data ...")
        print(NIc)
        print(Ic)
        message("Creating a design matrix ...")
        pd.desNIc <- NIc %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
          as.data.frame() %>% 
          dplyr::select(V8,V10, V11) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate( Genotype = gsub("[[:punct:]]*", "", Genotype),
                  Genotype = as.character(Genotype),
                  KP = ifelse(Genotype == "KP", 1, 0),
                  P = ifelse(Genotype == "P", 1, 0),
                  Time = gsub("T", "", Time),
                  Time = as.numeric(Time))
        pd.desIc <- Ic %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
          as.data.frame() %>% 
          dplyr::select(V8,V10, V11) %>%
          setNames(c("Genotype","Time", "Replicate"))%>%
          mutate( Genotype = gsub("[[:punct:]]*", "", Genotype),
                  Genotype = as.character(Genotype),
                  KP = ifelse(Genotype == "KP", 1, 0),
                  P = ifelse(Genotype == "P", 1, 0),
                  Time = gsub("T", "", Time),
                  Time = as.numeric(Time))
        pd.desc<- bind_rows(pd.desNIc,pd.desIc)
        pd.desc$Replicate<- rep(1:10, 1, each = 3)
        DT::datatable(head(pd.desc), caption = "Design matrix")
        rownames(pd.desc) <- pd$geo_accession[67:96]
        require(maSigPro)
        designc <- make.design.matrix(pd.desc[, c(2,3,4,5)], degree = 2)
        designc$edesign %>%
          as.data.frame() %>% 
          DT::datatable()
        message("Defining the regression model...")
        designc$groups.vector
        message("sourcing functions from github...")
        require(doParallel)
        require(itertools)
        require(devtools)
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/Tfitp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/pvectorp.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/position.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2sym.R")
        source_url("https://raw.githubusercontent.com/HediaTnani/Leishmania/master/probe2entrez.R")
        message("Fitting the regression model...")
        fitc <- p.vectorp(expc, designc,Q = 0.1, MT.adjust = "BH", min.obs = 3)
        tstepc <- T.fitp(fitc, step.method = "backward", alfa = 0.05)
        sigsc <- get.siggenes(tstepc, rsq = 0.6, vars = "groups")
        sigsc %>%
          glimpse()
        names(sigsc)
        names(sigsc$sig.genes)
        suma2Venn(sigsc$summary[,c(1:2)])
        clusters=see.genes(sigsc$sig.genes$PvsKP, show.fit = T, dis=designc$dis, cluster.method = "hclust", cluster.data = 1, k.mclust = TRUE)
        df=clusters$cut %>% as.data.frame() 
        write.csv(df, "PvsKP_C_Geneclusters.csv")
        PvsKP <- sigsc$sig.genes$PvsKP$sig.pvalues
        message ("from probes to entrez id ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                           dataset = "mmusculus_gene_ensembl",
                           host = 'www.ensembl.org')
        PvsKP <- probe2entrez(PvsKP, mart.mm)
        write.csv(PvsKP, "PvsKP_C_v0.csv")
        message ("from probes to genes symbols ....")
        require(biomaRt)
        mart.mm <- useMart(biomart = "ensembl",
                         dataset = "mmusculus_gene_ensembl",
                         host="https://dec2021.archive.ensembl.org")
        mart.hs = useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                        host="https://dec2021.archive.ensembl.org")
        PvsKP <- probe2sym(PvsKP, mart.mm, mart.hs)
        message ("writing csv files ....")
        write.csv(PvsKP, "PvsKP_C.csv")
        message ("Analysis PvsKP for c57 completed successfully ....")
      }#else
    }#PvsKP
  }#c57
}#function
