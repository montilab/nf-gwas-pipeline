#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
gds.file <- args[1]
nullmod <- args[2]
max.miss <- as.numeric(args[3])
result.file <- args[4]
log.file <- args[5]

sink(log.file, append=FALSE, split=TRUE)
date()
suppressPackageStartupMessages(library(GMMAT))

####Null model
nullmod <- readRDS(nullmod)

####Mixed effect model
cat("\n####assocTestAggregate starts\n")
glmm.score(nullmod, infile=gds.file, outfile=result.file)
cat("####assocTestAggregate ends\n\n")

date()
sink()
