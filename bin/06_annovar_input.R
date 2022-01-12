#!/usr/bin/env Rscript
args<-commandArgs(TRUE)
max_pval <- as.numeric(args[1])

sink('annovar_input.log', append=FALSE, split=TRUE)
date()
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))

snps <- fread("all_chr_caf_annotated.csv", header=T, stringsAsFactors=F)
snps <- as.tbl(snps)

snps <- snps %>%
  rename(pval = contains("pval"))

snps <- snps %>% filter(!is.na(pval))
print(paste0("number of variants analyzed : ",dim(snps)[1]))

snps <- snps %>% filter(pval < max_pval & pval != 0)
print(paste0("number of variants with p-value < ", max_pval," : ",dim(snps)[1]))

fwrite(snps, "top_snps_caf_annotated.csv", row.names = F)

snps.annovar <- snps %>%
  mutate(pos2 = pos, A1 = 0, A2 = 0) %>%
  select(chr, pos, pos2, A1, A2, snpID)

write.table(snps.annovar, "top_snps_input.txt", sep = "\t", row.names = F, col.names = F)

date()
sink()
