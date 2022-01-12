#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
n <- length(args)
gds.file <- args[n-6]
pheno.file <- args[n-5]
phenotypes <- args[n-4]
covariates <- unlist(strsplit(args[n-3], ","))
model <- args[n-2]
grm <- args[n-1]
slope <- args[n]

sink("nullmod_longitudinal.log", append=FALSE, split=TRUE)
date()
suppressPackageStartupMessages(library(SeqArray))
suppressPackageStartupMessages(library(GENESIS))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(SeqVarTools))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(SNPRelate))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GMMAT))

####Open GDS
gds <- seqOpen(gds.file)

####Create a SeqVarData object
pheno.dat <- read.csv(pheno.file, stringsAsFactors=F, header=T)
pheno.dat$sample.id <- as.character(pheno.dat[,1])
print(paste0("number of individuals in pheno data : ",length(pheno.dat$sample.id)))
if(length(pheno.dat$sample.id)>length(unique(pheno.dat$sample.id))){
    pheno.dat1 <- pheno.dat
    pheno.dat <- pheno.dat[!duplicated(pheno.dat$sample.id),]
}

gds.sample.id <- data.frame(sample.id=seqGetData(gds, "sample.id"),stringsAsFactors=F)
print(paste0("number of individuals in geno data : ", dim(gds.sample.id)[1]))

annot <- left_join(gds.sample.id, pheno.dat)

annot <- annot[,c("sample.id",colnames(annot)[colnames(annot)%in%c(phenotypes,covariates)])]

analysis.sample.id <- na.omit(annot[,c("sample.id",colnames(annot)[colnames(annot)%in%c(phenotypes,covariates)])])$sample.id
print(paste0("number of individuals to analyze : ",length(analysis.sample.id)))
saveRDS(analysis.sample.id, "analysis.sample.id.rds")

model.dat <- annot[annot$sample.id%in%analysis.sample.id,]

metadata <- data.frame(labelDescription=colnames(annot),row.names=names(annot))
annot <- AnnotatedDataFrame(annot, metadata)
saveRDS(annot, "annot.rds")

all.equal(annot$sample.id,seqGetData(gds, "sample.id"))
seqData <- SeqVarData(gds, sampleData=annot)

model.dat <- model.dat[, colnames(model.dat)%in%c("sample.id", paste0("PC",1:32))]
if(is.null(dim(model.dat))){
    pheno.pc.dat <- pheno.dat1[pheno.dat1$sample.id%in%analysis.sample.id,]
}else{
    pheno.dat1 <- pheno.dat1[pheno.dat1$sample.id%in%analysis.sample.id,]
    pheno.pc.dat <- left_join(pheno.dat1, model.dat, by="sample.id")
}

####Mixed effect null model
if(model=="linear"){
	model.switch <- gaussian(link = "identity")
}
if(model=="logistic"){
	model.switch <- binomial(link = "logit")
}

fix.eff=paste(phenotypes,"~ 1")
if(!is.null(covariates)){
	for(covi in covariates)
		fix.eff=paste(fix.eff,"+",covi)
}
fix.eff=formula(fix.eff)

if(grm=="null"){
	if(slope=="null"){
		cat("\n####glmmkin starts\n")
		nullmod <- glmmkin(fix.eff, data=pheno.pc.dat, kins=NULL, 
	        		   id="sample.id", family = model.switch)
		cat("####glmmkin ends\n\n")
	}else{
		cat("\n####glmmkin with grm starts\n")
		nullmod <- glmmkin(fix.eff, data=pheno.pc.dat, kins=NULL, 
		           random.slope=slope, id="sample.id", family = model.switch)
		cat("####glmmkin with grm ends\n\n")
	}	
}else{
	grm <- readRDS(grm)
		if(slope=="null"){
		cat("\n####glmmkin starts\n")
		nullmod <- glmmkin(fix.eff, data=pheno.pc.dat, kins=as.matrix(grm), 
	        		   id="sample.id", family = model.switch)
		cat("####glmmkin ends\n\n")
	}else{
		cat("\n####glmmkin starts\n")
		nullmod <- glmmkin(fix.eff, data=pheno.pc.dat, kins=as.matrix(grm), 
		           random.slope=slope, id="sample.id", family = model.switch)
		cat("####glmmkin ends\n\n")
	}
}

saveRDS(nullmod,"nullmod_longitudinal.rds")

date()
sink()
