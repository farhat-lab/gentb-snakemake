#!/usr/bin/env Rscript

## Minimal prediction routine with all possible preprocessing
## Configured for calling by Rscript
## Call as:
## Rscript TBpredict.R "testinputfile.csv"

arg <- commandArgs(trailingOnly = TRUE)

library(foreign)
library(jsonlite)
library(randomForest)

predictfunction <- function(filename){
set.seed(5414)
#input the variant matrix
strain <- read.csv('test_pyra_script.matrix', header=TRUE)

#prepare output for variants
important_strain<-matrix(NA,nrow=5,ncol=1,dimnames=list(c(), 'pza'))
important<-vector("list", nrow(strain))
names(important)<-strain[,1]

#prepare output for RandomForest probability
result <- matrix(NA, nrow=1, ncol=3)

#load RandomForest object
load("pza_finalpredict.RData")

#perform prediction and write strain ID, drug, and probability to result
Valid <- predict(drugg.full.rf, strain, type='prob', norm.votes=TRUE, predict.all=FALSE)
result[1,] <- c(as.vector(strain[1,1]), 'pza', round(Valid[1,1],3))

#write the variant used for prediction to object
imp<-colnames(strain)[which(strain[1,]==1)]
important_strain[1:length(imp),1]<-imp
important[[1]]<-important_strain

#bind RandomForest probability and variant used to final result object
l <- list(result, important)

## Save JSON file
file_noext <- substr(filename, 1, nchar(filename) - 4)
cat(toJSON(l, pretty = TRUE), "\n", file = paste0(file_noext, ".json"))
}

suppressWarnings(predictfunction(arg))