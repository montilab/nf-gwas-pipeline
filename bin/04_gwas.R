#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
gds.file <- args[1]
annot <- args[2]
nullmod <- args[3]
test <- args[4]
imputed <- args[5]
result.file <- args[6]
log.file <- args[7]

sink(log.file, append=FALSE, split=TRUE)
date()
suppressPackageStartupMessages(library(SeqArray))
suppressPackageStartupMessages(library(GENESIS))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(SeqVarTools))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))

####Open GDS
gds <- seqOpen(gds.file)

#save vector with snps rs ids
snps <- data.frame(variant.id=seqGetData(gds, "variant.id"), snpID=seqGetData(gds, "annotation/id"))

####Create a SeqVarData object
annot <- readRDS(annot)
seqData <- SeqVarData(gds, sampleData=annot)

####Null model
nullmod <- readRDS(nullmod)

####GWAS
iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)

cat("\n####assocTestSingle starts\n")
assoc <- assocTestSingle(iterator, nullmod, test=test, imputed=ifelse(imputed=="true",T,F), verbose=T)
cat("####assocTestSingle ends\n\n")

assoc <- left_join(assoc, snps, by="variant.id")
colnames(assoc)[colnames(assoc)=="n.obs"] <- "N"

fwrite(assoc, file = result.file, row.names=FALSE)

date()
sink()
