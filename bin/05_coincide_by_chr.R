#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
in.file <- args[1]
out.file <- args[2] 
log.file <- args[3]

sink(log.file, append=FALSE, split=TRUE)
date()
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))

result.dat <- fread(in.file, sep="\t", stringsAsFactors = F, header = T)

colnames(result.dat)[colnames(result.dat)=="SNP"] <- "snpID"
colnames(result.dat)[colnames(result.dat)=="CHR"] <- "chr"
colnames(result.dat)[colnames(result.dat)=="POS"] <- "pos"
colnames(result.dat)[colnames(result.dat)=="REF"] <- "ref"
colnames(result.dat)[colnames(result.dat)=="ALT"] <- "alt"
colnames(result.dat)[colnames(result.dat)=="SCORE"] <- "Score"
colnames(result.dat)[colnames(result.dat)=="PVAL"] <- "Score.pval"

result.dat$Score.pval <- as.numeric(result.dat$Score.pval)
result.dat$Score.SE <- sqrt(result.dat$VAR)

fwrite(result.dat, out.file, row.names=F)

date()
sink()
