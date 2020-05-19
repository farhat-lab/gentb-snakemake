#!/usr/bin/env Rscript

## Script to run as onstart script to install the required R packages.

local({r<-getOption("repos")
      r["CRAN"]<-"http://cran.fhcrc.org/"
      options(repos=r)})

script.name <- sub(prefix, "", options[grep(prefix, options)])
script.basename <- dirname(script.name)
libs <- file.path(script.basename, '../data/Rlib')

install.packages("jsonlite", quiet=TRUE, lib=libs)
install.packages("randomForest", quiet=TRUE, lib=libs)