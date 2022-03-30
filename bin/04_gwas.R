#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
gds.file <- args[1]
annot1 <- args[2]
annot2 <- args[3]
nullmod <- args[4]
max.miss <- as.numeric(args[5])
test <- args[6]
imputed <- args[7]
result.file <- args[8]
log.file <- args[9]

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
try(annot <- readRDS(annot1))
try(annot <- readRDS(annot2))
seqData <- SeqVarData(gds, sampleData=annot)

####Null model
nullmod <- readRDS(nullmod)

####Calculate Missing Rate of SNPs
seqSetFilter(seqData, sample.id=nullmod$fit$sample.id)
snps$missing.rate <- seqMissing(seqData)

snps.clean <- snps[snps$missing.rate < (1-max.miss), ]

####GWAS
seqSetFilter(seqData, variant.id=snps.clean$variant.id)

iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)

cat("\n####assocTestSingle starts\n")
assoc <- assocTestSingle(iterator, nullmod, test=test, imputed=ifelse(imputed=="true",T,F), verbose=T)
cat("####assocTestSingle ends\n\n")

assoc <- left_join(assoc, snps.clean, by="variant.id")
colnames(assoc)[colnames(assoc)=="n.obs"] <- "N"

fwrite(assoc, file = result.file, row.names=FALSE)

date()
sink()
