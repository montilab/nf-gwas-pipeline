#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
pheno.file <- args[1]
phenotypes <- args[2]
covariates <- unlist(strsplit(args[3], ","))
model <- args[4]
analysis.sample.id <- args[5]
pc_df <- args[6]
grm <- args[7]
slope <- args[8]

sink("nullmod_longitudinal.log", append=FALSE, split=TRUE)
date()
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GMMAT))
suppressPackageStartupMessages(library(dplyr))

####Phenotype data
pc.df <- readRDS(pc_df)
analysis.sample.id <- readRDS(analysis.sample.id)
pc.df <- pc.df[pc.df$sample.id%in%analysis.sample.id,]
pc.df <- pc.df[, colnames(pc.df)%in%c("sample.id", paste0("PC",1:32))]
pheno.dat <- read.csv(pheno.file, stringsAsFactors=F, header=T)
pheno.dat$sample.id <- as.character(pheno.dat[,1])

if(is.null(dim(pc.df))){
    pheno.pc.dat <- pheno.dat[pheno.dat$sample.id%in%analysis.sample.id,]
}else{
    pheno.dat <- pheno.dat[pheno.dat$sample.id%in%analysis.sample.id,]
    pheno.pc.dat <- left_join(pheno.dat, pc.df, by="sample.id")
}

####covariance matrix from pcrelate output
grm <- readRDS(grm)

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

if(slope=="null"){
	cat("\n####glmmkin starts\n")
	nullmod <- glmmkin(fix.eff, data=pheno.pc.dat, kins=as.matrix(grm),
	           random.slope=NULL, id="sample.id", family = model.switch)
	cat("####glmmkin ends\n\n")

}else{
	cat("\n####glmmkin starts\n")
	nullmod <- glmmkin(fix.eff, data=pheno.pc.dat, kins=as.matrix(grm), 
	           random.slope=slope, id="sample.id", family = model.switch)
	cat("####glmmkin ends\n\n")

}

saveRDS(nullmod,"nullmod_longitudinal.rds")

date()
sink()
