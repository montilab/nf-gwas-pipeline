#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

sink('combine_results.log', append=FALSE, split=TRUE)
date()
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))

result <- fread(args[1], stringsAsFactors=F, header=T)
for (i in 2:length(args)) {
  res <- fread(args[i], stringsAsFactors=F, header=T)
  result <- rbind(result, res)
}

fwrite(result, paste0("all_chr_caf_annotated.csv"), row.names = FALSE)

date()
sink()
