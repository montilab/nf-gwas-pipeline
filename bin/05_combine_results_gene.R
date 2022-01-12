#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path1 <- args[1]

sink('combine_results.log', append=FALSE, split=TRUE)
date()
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))

result <- fread(args[1], stringsAsFactors=F, header=T)
try(for (i in 2:length(args)) {
      res <- fread(args[i], stringsAsFactors=F, header=T)
      result <- rbind(result, res)
    })

result <- result %>%
  rename(pval = contains("pval"))

result <- result %>%
	select ("chr", "start", "end", "width", "strand", "gene_id", "n.site", "n.alt", "n.sample.alt", contains("Score"), contains("Wald"), "pval")

fwrite(result[order(result$pval),], paste0("all_chr.csv"), row.names = FALSE)

date()
sink()
