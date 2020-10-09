#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
ref <- args[1]

sink('add_anno_results.log', append=FALSE, split=TRUE)
date()
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))

results <- fread("top_snps_caf_annotated.csv", header=T, stringsAsFactors=F)
annovar <- fread(paste0("top_annotation.",ref,"_multianno.csv"), header=T, stringsAsFactors=F)

annovar$chr <- annovar$Chr
annovar$pos <- annovar$Start

annot.results <- left_join(results, annovar, by = c("chr", "pos") )

annot.results <- annot.results %>%
	select ("snpID", "chr", "pos", "REF", "ALT", contains("Imputation_Rsq"), contains("Imputation_mark"), contains("Score"), contains("Wald"), "pval", "N", contains("n.case"), contains("n.control"), contains("caf"), contains("dosage"), contains("refGene"))

fwrite(annot.results[order(annot.results$pval),], "top_snps_annotation.csv", row.names = FALSE)

date()
sink()
