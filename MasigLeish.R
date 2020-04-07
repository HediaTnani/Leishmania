# Author: Hédia Tnani
# Date: 3/4/2020
# Usage:MasigLeish(mouse = "Balb") or MasigLeish(mouse = "c57")
##############################################################################################################################
MasigLeish = function (mouse){ 
  options(warn = -1)
  cat ("Install required packages")
  packages <- c("tidyverse", "doParallel","foreign","itertools", "devtools")
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
  if (mouse== "Balb"){# Balb
    message ("loading packages ....")
    require(GEOquery)
    require(oligo)
    require(mogene10sttranscriptcluster.db)
    require(RCurl)
    require(foreign)
    require(tidyverse)
    setwd(getwd())
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
    message("Annotation ...")
    annotation(eset_rma) <- "mogene10sttranscriptcluster.db"
    keytypes(mogene10sttranscriptcluster.db)
    gns <- select(mogene10sttranscriptcluster.db, keys(mogene10sttranscriptcluster.db), c("ENTREZID", "SYMBOL"))
    gns <- gns[!duplicated(gns[,1]),]
    gns = gns[,-1]
    row.names(gns) = keys(mogene10sttranscriptcluster.db)
    #exp <- exprs(eset_rma)
    exp <- data.frame(exprs(eset_rma))
    exp.anno <- merge(gns, exp, by.x=0, by.y=0, all=TRUE)
    expb <- exp[,1:48]
    message("Viewing expression data ...")
    expb %>%
      as.data.frame() %>%
      head()%>%
      dplyr::select(1:4) 
    message("Getting phenotypic data ...")
    gse <- getGEO("GSE31997", AnnotGPL = TRUE)
    pd <- pData(gse$GSE31997_series_matrix.txt.gz)
    NIb= as.character(pd$title[1:18]) 
    Ib= as.character(pd$title[19:48]) 
    message("Printing phenotypic data ...")
    print(NIb)
    print(Ib)
    message("Creating a design matrix ...")
    pd.desNIb <- NIb %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
      as.data.frame() %>% 
      unite("V2_V3", c("V2","V3"),sep = " " ) %>%
      dplyr::select(V2_V3, V4, V7, V8) %>%
      setNames(c("Genotype", "Type", "Time", "Replicate"))%>%
      mutate(Control = ifelse(Genotype == "Non Infected", 1, 0),
               Genotype = gsub("[[:punct:]]*", "", Genotype),
               Genotype = as.character(Genotype),
               P = ifelse(Genotype == "P", 1, 0),
               KP = ifelse(Genotype == "KP", 1, 0),
               Time = gsub("T", "", Time),
               Time = as.numeric(Time))
    pd.desIb <- Ib %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
      as.data.frame() %>% 
      dplyr::select(V8, V2, V10, V11) %>%
      setNames(c("Genotype", "Type", "Time", "Replicate"))%>%
      mutate( Control = ifelse(Genotype == "Non Infected", 1, 0),
               Genotype = gsub("[[:punct:]]*", "", Genotype),
               Genotype = as.character(Genotype),
               P = ifelse(Genotype == "P", 1, 0),
               KP = ifelse(Genotype == "KP", 1, 0),
               Time = gsub("T", "", Time),
               Time = as.numeric(Time))
    pd.desb<- bind_rows(pd.desNIb,pd.desIb)
    pd.desb$Replicate<- rep(1:16, 1, each = 3)
    DT::datatable(head(pd.desb), caption = "Design matrix")
    rownames(pd.desb) <- pd$geo_accession[1:48]
    require(maSigPro)
    designb <- make.design.matrix(pd.desb[, c(3,4,5,6,7)], degree = 2)
    designb$edesign %>%
      as.data.frame() %>% 
      DT::datatable()
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
    message("Fitting the regression model...")
    fitb<- p.vectorp(expb, designb,Q = 0.05, MT.adjust = "BH", min.obs = 3)
    tstepb <- T.fitp(fitb, step.method = "backward", alfa = 0.05)
    sigsb <- get.siggenes(tstepb, rsq = 0.6, vars = "groups")
    sigsb %>%
      glimpse()
    control   <- sigsb$sig.genes$Control$sig.pvalues
    P <- sigsb$sig.genes$PvsControl$sig.pvalues
    KP <- sigsb$sig.genes$KPvsControl$sig.pvalues
    message ("from probes to genes symbols ....")
    require(biomaRt)
    mart.mm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                       dataset = "mmusculus_gene_ensembl", 
                       host = 'www.ensembl.org')
    control   <- probe2sym(control, mart.mm)
    P <- probe2sym(P, mart.mm)
    KP <- probe2sym(KP, mart.mm)
    message ("writing csv files ....")
    write.csv(control, "control_Balb.csv")
    write.csv(P, "P_Balb.csv")
    write.csv(KP, "KP_Balb.csv")
    return (control)
    return(P)
    return (KP)
    message ("Analysis completed successfully ....")
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
      message("Annotation ...")
      annotation(eset_rma) <- "mogene10sttranscriptcluster.db"
      keytypes(mogene10sttranscriptcluster.db)
      gns <- select(mogene10sttranscriptcluster.db, keys(mogene10sttranscriptcluster.db), c("ENTREZID", "SYMBOL"))
      gns <- gns[!duplicated(gns[,1]),]
      gns = gns[,-1]
      row.names(gns) = keys(mogene10sttranscriptcluster.db)
      #exp <- exprs(eset_rma)
      exp <- data.frame(exprs(eset_rma))
      exp.anno <- merge(gns, exp, by.x=0, by.y=0, all=TRUE)
      expb <- exp[,1:48]
      message("Viewing expression data ...")
      expb %>%
        as.data.frame() %>%
        head()%>%
        dplyr::select(1:4) 
      message("Getting phenotypic data ...")
      gse <- getGEO("GSE31997", AnnotGPL = TRUE)
      pd <- pData(gse$GSE31997_series_matrix.txt.gz)
      NIb= as.character(pd$title[1:18]) 
      Ib= as.character(pd$title[19:48]) 
      message("Printing phenotypic data ...")
      print(NIb)
      print(Ib)
      message("Creating a design matrix ...")
      pd.desNIb <- NIb %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
        as.data.frame() %>% 
        unite("V2_V3", c("V2","V3"),sep = " " ) %>%
        dplyr::select(V2_V3, V4, V7, V8) %>%
        setNames(c("Genotype", "Type", "Time", "Replicate"))%>%
        mutate(Control = ifelse(Genotype == "Non Infected", 1, 0),
               Genotype = gsub("[[:punct:]]*", "", Genotype),
               Genotype = as.character(Genotype),
               P = ifelse(Genotype == "P", 1, 0),
               KP = ifelse(Genotype == "KP", 1, 0),
               Time = gsub("T", "", Time),
               Time = as.numeric(Time))
      pd.desIb <- Ib %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
        as.data.frame() %>% 
        dplyr::select(V8, V2, V10, V11) %>%
        setNames(c("Genotype", "Type", "Time", "Replicate"))%>%
        mutate( Control = ifelse(Genotype == "Non Infected", 1, 0),
                Genotype = gsub("[[:punct:]]*", "", Genotype),
                Genotype = as.character(Genotype),
                P = ifelse(Genotype == "P", 1, 0),
                KP = ifelse(Genotype == "KP", 1, 0),
                Time = gsub("T", "", Time),
                Time = as.numeric(Time))
      pd.desb<- bind_rows(pd.desNIb,pd.desIb)
      pd.desb$Replicate<- rep(1:16, 1, each = 3)
      DT::datatable(head(pd.desb), caption = "Design matrix")
      rownames(pd.desb) <- pd$geo_accession[1:48]
      require(maSigPro)
      designb <- make.design.matrix(pd.desb[, c(3,4,5,6,7)], degree = 2)
      designb$edesign %>%
        as.data.frame() %>% 
        DT::datatable()
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
      message("Fitting the regression model...")
      fitb<- p.vectorp(expb, designb,Q = 0.05, MT.adjust = "BH", min.obs = 3)
      tstepb <- T.fitp(fitb, step.method = "backward", alfa = 0.05)
      sigsb <- get.siggenes(tstepb, rsq = 0.6, vars = "groups")
      sigsb %>%
        glimpse()
      control   <- sigsb$sig.genes$Control$sig.pvalues
      P <- sigsb$sig.genes$PvsControl$sig.pvalues
      KP <- sigsb$sig.genes$KPvsControl$sig.pvalues
      message ("from probes to genes symbols ....")
      require(biomaRt)
      mart.mm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                         dataset = "mmusculus_gene_ensembl", 
                         host = 'www.ensembl.org')
      control   <- probe2sym(control, mart.mm)
      P <- probe2sym(P, mart.mm)
      KP <- probe2sym(KP, mart.mm)
      message ("writing csv files ....")
      write.csv(control, "control_Balb.csv")
      write.csv(P, "P_Balb.csv")
      write.csv(KP, "KP_Balb.csv")
      return (control)
      return(P)
      return (KP)
      message ("Analysis completed successfully ....")
      }# else
      }#Balb
  
  if(mouse=="c57"){#c57
    message (" loading packages ....")
    require(GEOquery)
    require(oligo)
    require(mogene10sttranscriptcluster.db)
    require(RCurl)
    require(foreign)
    require(tidyverse)
    setwd(getwd())
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
    message("Annotation ...")
    annotation(eset_rma) <- "mogene10sttranscriptcluster.db"
    keytypes(mogene10sttranscriptcluster.db)
    gns <- select(mogene10sttranscriptcluster.db, keys(mogene10sttranscriptcluster.db), c("ENTREZID", "SYMBOL"))
    gns <- gns[!duplicated(gns[,1]),]
    gns = gns[,-1]
    row.names(gns) = keys(mogene10sttranscriptcluster.db)
    #exp <- exprs(eset_rma)
    exp <- data.frame(exprs(eset_rma))
    exp.anno <- merge(gns, exp, by.x=0, by.y=0, all=TRUE)
    expc <- exp[,49:96]
    message("Viewing expression data ...")
    expc %>%
      as.data.frame() %>%
      head()%>%
      dplyr::select(1:4) 
    message("Getting phenotypic data ...")
    gse <- getGEO("GSE31997", AnnotGPL = TRUE)
    pd <- pData(gse$GSE31997_series_matrix.txt.gz)
    NIc= as.character(pd$title[49:66]) 
    Ic= as.character(pd$title[67:96]) 
    message("Printing phenotypic data ...")
    print(NIc)
    print(Ic)
    message("Creating a design matrix ...")
    pd.desNIc <- NIc %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
      as.data.frame() %>% 
      unite("V2_V3", c("V2","V3"),sep = " " ) %>%
      dplyr::select(V2_V3, V4, V7, V8) %>%
      setNames(c("Genotype", "Type", "Time", "Replicate"))%>%
      mutate(  Control = ifelse(Genotype == "Non Infected", 1, 0),
               Genotype = gsub("[[:punct:]]*", "", Genotype),
               Genotype = as.character(Genotype),
               P = ifelse(Genotype == "P", 1, 0),
               KP = ifelse(Genotype == "KP", 1, 0),
               Time = gsub("T", "", Time),
               Time = as.numeric(Time))
    
    pd.desIc <- Ic %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
      as.data.frame() %>% 
      dplyr::select(V8, V2, V10, V11) %>%
      setNames(c("Genotype", "Type", "Time", "Replicate"))%>%
      mutate(  Control = ifelse(Genotype == "Non Infected", 1, 0),
               Genotype = gsub("[[:punct:]]*", "", Genotype),
               Genotype = as.character(Genotype),
               P = ifelse(Genotype == "P", 1, 0),
               KP = ifelse(Genotype == "KP", 1, 0),
               Time = gsub("T", "", Time),
               Time = as.numeric(Time))
    pd.desc<- bind_rows(pd.desNIc,pd.desIc)
    pd.desc$Replicate<- rep(1:16, 1, each = 3)
    DT::datatable(head(pd.desc), caption = "Design matrix")
    rownames(pd.desc) <- pd$geo_accession[49:96]
    require(maSigPro)
    designc <- make.design.matrix(pd.desc[, c(3,4,5,6,7)], degree = 2)
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
    message("Fitting the regression model...")
    fitc<- p.vectorp(expc, designc,Q = 0.05, MT.adjust = "BH", min.obs = 3)
    tstepc <- T.fitp(fitc, step.method = "backward", alfa = 0.05)
    sigsc <- get.siggenes(tstepc, rsq = 0.6, vars = "groups")
    sigsc %>%
      glimpse()
    control   <- sigsc$sig.genes$Control$sig.pvalues
    P <- sigsc$sig.genes$PvsControl$sig.pvalues
    KP <- sigsc$sig.genes$KPvsControl$sig.pvalues
    message ("from probes to genes symbols ....")
    require(biomaRt)
    mart.mm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                       dataset = "mmusculus_gene_ensembl", 
                       host = 'www.ensembl.org')
    control   <- probe2sym(control, mart.mm)
    P <- probe2sym(P, mart.mm)
    KP <- probe2sym(KP, mart.mm)
    message ("writing csv files ....")
    write.csv(control, "control_C57.csv")
    write.csv(P, "P_C57.csv")
    write.csv(KP, "KP_C57.csv")
    return (control)
    return(P)
    return (KP)
    message ("Analysis completed successfully ....")
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
      message("Annotation ...")
      annotation(eset_rma) <- "mogene10sttranscriptcluster.db"
      keytypes(mogene10sttranscriptcluster.db)
      gns <- select(mogene10sttranscriptcluster.db, keys(mogene10sttranscriptcluster.db), c("ENTREZID", "SYMBOL"))
      gns <- gns[!duplicated(gns[,1]),]
      gns = gns[,-1]
      row.names(gns) = keys(mogene10sttranscriptcluster.db)
      #exp <- exprs(eset_rma)
      exp <- data.frame(exprs(eset_rma))
      exp.anno <- merge(gns, exp, by.x=0, by.y=0, all=TRUE)
      expc <- exp[,49:96]
      message("Viewing expression data ...")
      expc %>%
        as.data.frame() %>%
        head()%>%
        dplyr::select(1:4) 
      message("Getting phenotypic data ...")
      gse <- getGEO("GSE31997", AnnotGPL = TRUE)
      pd <- pData(gse$GSE31997_series_matrix.txt.gz)
      NIc= as.character(pd$title[49:66]) 
      Ic= as.character(pd$title[67:96]) 
      message("Printing phenotypic data ...")
      print(NIc)
      print(Ic)
      message("Creating a design matrix ...")
      pd.desNIc <- NIc %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
        as.data.frame() %>% 
        unite("V2_V3", c("V2","V3"),sep = " " ) %>%
        dplyr::select(V2_V3, V4, V7, V8) %>%
        setNames(c("Genotype", "Type", "Time", "Replicate"))%>%
        mutate(  Control = ifelse(Genotype == "Non Infected", 1, 0),
                 Genotype = gsub("[[:punct:]]*", "", Genotype),
                 Genotype = as.character(Genotype),
                 P = ifelse(Genotype == "P", 1, 0),
                 KP = ifelse(Genotype == "KP", 1, 0),
                 Time = gsub("T", "", Time),
                 Time = as.numeric(Time))
      
      pd.desIc <- Ic %>% str_match("^(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?)\\s+(.*?),\\s+(.*?)$") %>% 
        as.data.frame() %>% 
        dplyr::select(V8, V2, V10, V11) %>%
        setNames(c("Genotype", "Type", "Time", "Replicate"))%>%
        mutate(  Control = ifelse(Genotype == "Non Infected", 1, 0),
                 Genotype = gsub("[[:punct:]]*", "", Genotype),
                 Genotype = as.character(Genotype),
                 P = ifelse(Genotype == "P", 1, 0),
                 KP = ifelse(Genotype == "KP", 1, 0),
                 Time = gsub("T", "", Time),
                 Time = as.numeric(Time))
      pd.desc<- bind_rows(pd.desNIc,pd.desIc)
      pd.desc$Replicate<- rep(1:16, 1, each = 3)
      DT::datatable(head(pd.desc), caption = "Design matrix")
      rownames(pd.desc) <- pd$geo_accession[49:96]
      require(maSigPro)
      designc <- make.design.matrix(pd.desc[, c(3,4,5,6,7)], degree = 2)
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
      message("Fitting the regression model...")
      fitc<- p.vectorp(expc, designc,Q = 0.05, MT.adjust = "BH", min.obs = 3)
      tstepc <- T.fitp(fitc, step.method = "backward", alfa = 0.05)
      sigsc <- get.siggenes(tstepc, rsq = 0.6, vars = "groups")
      sigsc %>%
        glimpse()
      control   <- sigsc$sig.genes$Control$sig.pvalues
      P <- sigsc$sig.genes$PvsControl$sig.pvalues
      KP <- sigsc$sig.genes$KPvsControl$sig.pvalues
      message ("from probes to genes symbols ....")
      require(biomaRt)
      mart.mm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                         dataset = "mmusculus_gene_ensembl", 
                         host = 'www.ensembl.org')
      control   <- probe2sym(control, mart.mm)
      P <- probe2sym(P, mart.mm)
      KP <- probe2sym(KP, mart.mm)
      message ("writing csv files ....")
      write.csv(control, "control_C57.csv")
      write.csv(P, "P_C57.csv")
      write.csv(KP, "KP_C57.csv")
      return (control)
      return(P)
      return (KP)
      message ("Analysis completed successfully ....")
    }# else
    }#c57
  
}# function

