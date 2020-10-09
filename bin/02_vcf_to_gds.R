#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

vcf.file <- args[1]
gds.file <- args[2]
log.file <- args[3]

sink(log.file, append=FALSE, split=TRUE)
date()
suppressPackageStartupMessages(library(SeqArray))

cat("\n####seqVCF2GDS starts\n")
seqVCF2GDS(vcf.file, gds.file, verbose=T)
cat("####seqVCF2GDS ends\n\n")

date()
sink()
