#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

sink("merge_gds.log", append=FALSE, split=TRUE)
date()
suppressPackageStartupMessages(library(SeqArray))

cat("\n####seqMerge starts\n")
seqMerge(sort(args), "merged.gds")
cat("####seqMerge ends\n\n")

date()
sink()
