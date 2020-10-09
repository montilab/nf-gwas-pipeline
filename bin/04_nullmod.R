#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
gds.file <- args[1]
phenotypes <- args[2]
covariates <- unlist(strsplit(args[3], ","))
model <- args[4]
analysis.sample.id <- args[5]
annot <- args[6]
pc_df <- args[7]
grm <- args[8]

sink("nullmod.log", append=FALSE, split=TRUE)
date()
suppressPackageStartupMessages(library(SeqArray))
suppressPackageStartupMessages(library(GENESIS))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(SeqVarTools))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(SNPRelate))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))

####Open GDS
gds <- seqOpen(gds.file)

####Create a SeqVarData object
analysis.sample.id <- readRDS(analysis.sample.id)
annot <- readRDS(annot)
seqData <- SeqVarData(gds, sampleData=annot)

####PCA
pc.df <- readRDS(pc_df)

####add PCs to sample annotation in SeqVarData object
annot <- AnnotatedDataFrame(pc.df)
saveRDS(annot, "annot_pc.rds")
sampleData(seqData) <- annot

####covariance matrix from pcrelate output
grm <- readRDS(grm)

####Null model
if(model=="linear"){
	model.switch <- "gaussian"
}
if(model=="logistic"){
	model.switch <- "binomial"
}

cat("\n####fitNullModel starts\n")
nullmod <- fitNullModel(seqData, outcome=phenotypes, 
                        covars=covariates,
                        cov.mat=grm,
                        family=model.switch, verbose=T)
cat("####fitNullModel ends\n\n")

saveRDS(nullmod,"nullmod.rds")

date()
sink()
