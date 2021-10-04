#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
n <- length(args)
gds.file <- args[n-5]
pheno.file <- args[n-4]
phenotypes <- args[n-3]
covariates <- unlist(strsplit(args[n-2], ","))
model <- args[n-1]
grm <- args[n]

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
pheno.dat <- read.csv(pheno.file, stringsAsFactors=F, header=T)
pheno.dat$sample.id <- as.character(pheno.dat[,1])
print(paste0("number of individuals in pheno data : ",length(pheno.dat$sample.id)))

gds.sample.id <- data.frame(sample.id=seqGetData(gds, "sample.id"),stringsAsFactors=F)
print(paste0("number of individuals in geno data : ", dim(gds.sample.id)[1]))

annot <- left_join(gds.sample.id, pheno.dat)

annot <- annot[,c("sample.id",colnames(annot)[colnames(annot)%in%c(phenotypes,covariates)])]

analysis.sample.id <- na.omit(annot[,c("sample.id",colnames(annot)[colnames(annot)%in%c(phenotypes,covariates)])])$sample.id
print(paste0("number of individuals to analyze : ",length(analysis.sample.id)))
saveRDS(analysis.sample.id, "analysis.sample.id.rds")

metadata <- data.frame(labelDescription=colnames(annot),row.names=names(annot))
annot <- AnnotatedDataFrame(annot, metadata)
saveRDS(annot, "annot.rds")

all.equal(annot$sample.id,seqGetData(gds, "sample.id"))
seqData <- SeqVarData(gds, sampleData=annot)

####Null model
if(model=="linear"){
	model.switch <- "gaussian"
}
if(model=="logistic"){
	model.switch <- "binomial"
}

if(grm=="null"){
	cat("\n####fitNullModel starts\n")
	nullmod <- fitNullModel(seqData, outcome=phenotypes, 
    	                    covars=covariates,
            	            family=model.switch, verbose=T)
	cat("####fitNullModel ends\n\n")
}else{
	grm <- readRDS(grm)
	cat("\n####fitNullModel with grm starts\n")
	nullmod <- fitNullModel(seqData, outcome=phenotypes, 
    	                    covars=covariates,
        	                cov.mat=grm,
            	            family=model.switch, verbose=T)
	cat("####fitNullModel with grm ends\n\n")
}

saveRDS(nullmod,"nullmod.rds")

date()
sink()
