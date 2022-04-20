#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
mac = as.numeric(args[1])

sink('combine_results.log', append=FALSE, split=TRUE)
date()
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))

result <- fread(args[2], stringsAsFactors=F, header=T)
try(result <- result %>% filter(MAC > mac))
try(for (i in 3:length(args)) {
      res <- fread(args[i], stringsAsFactors=F, header=T)
      try(res <- res %>% filter(MAC > mac))
      result <- rbind(result, res)
    })

fwrite(result, paste0("all_chr_caf_annotated.csv"), row.names = FALSE)

date()
sink()
