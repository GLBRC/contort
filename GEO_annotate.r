#! /Library/Frameworks/R.framework/Versions/Current/Resources/Rscript
# Parsing GEO Files to include normalized data and associated gene names

### Check for required packages and install if needed and then load libraries

CRAN_packages <- function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x, dependencies = TRUE, repos="http://cran.r-project.org", quiet = TRUE)
    require(x,character.only=TRUE)
  }
}

CRAN_packages(tidyverse)
CRAN_packages(BiocManager)

BIOCONDUCTOR_packages<-function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    BiocManager::install(x, ask = FALSE, quiet = TRUE)
    require(x,character.only=TRUE)
  }
}

BIOCONDUCTOR_packages(GEOquery)

library(GEOquery)
library(tidyverse)

### Set path and load GEOpatch.r IF R VERSION < 3.4

version = strsplit(sessionInfo()$R.version$version.string, split = " ")[[1]][3]

if (version < 3.4) {

  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  script.dirname <- dirname(script.name)
  sourcename = paste(script.dirname, "GEOpatch.r", sep = "/")

  source(sourcename)
}

### Setup
rm(list = ls())
setwd(".")
cwd = getwd()

# generate list of GEO files to use
gse_names <- read.table(file = "series.txt", header = FALSE)

# download GSE and GPL and combine all annotations with data files
for(i in 1:nrow(gse_names)){
  gse_name <- gse_names[i, 1]
  gse <- try(getGEO(gse_name, GSEMatrix = TRUE))
  if(inherits(gse, "try-error")) next
  data <- exprs(gse[[1]])
  pData <- pData(featureData(gse[[1]]))
  if (nrow(pData)>0){
    pData <- sapply(pData, function(f) {is.na(f) <- which(f == ""); f})
    combo <- cbind(pData, data)
    fileName <- paste(gse_name, "_GEOquery_results.txt")
    fileName <- gsub(" ", "", fileName, fixed = TRUE)
    write.table(combo, file = fileName, quote = FALSE, sep = "\t", row.names = FALSE)
  }
}
