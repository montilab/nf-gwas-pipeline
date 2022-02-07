#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
gds.file <- args[1]
pheno.file <- args[2]
phenotypes <- args[3]
covariates <- unlist(strsplit(args[4], ","))
snpset.file <- args[5]

sink(log.file, append=FALSE, split=TRUE)

date()
suppressPackageStartupMessages(library(SeqArray))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(SNPRelate))
suppressPackageStartupMessages(library(data.table))

####Open GDS
gds <- seqOpen(gds.file)

####Analysis sample set
pheno.dat <- read.csv(pheno.file, stringsAsFactors=F, header=T)
pheno.dat$sample.id <- as.character(pheno.dat[,1])
print(paste0("number of individuals in pheno data : ",length(pheno.dat$sample.id)))
if(length(pheno.dat$sample.id)>length(unique(pheno.dat$sample.id))){
  pheno.dat <- pheno.dat[!duplicated(pheno.dat$sample.id),]
}

for(i in 1:32){
  if(sum(colnames(pheno.dat)==paste0("PC",i))==1){
    colnames(pheno.dat)[colnames(pheno.dat)==paste0("PC",i)] <- paste0("PC",i,".pheno")
  }
}

gds.sample.id <- data.frame(sample.id=seqGetData(gds, "sample.id"),stringsAsFactors=F)
print(paste0("number of individuals in geno data : ", dim(gds.sample.id)[1]))

annot <- left_join(gds.sample.id, pheno.dat)

annot <- annot[,c("sample.id",colnames(annot)[colnames(annot)%in%c(phenotypes,covariates)])]

analysis.sample.id <- na.omit(annot[,c("sample.id",colnames(annot)[colnames(annot)%in%c(phenotypes,covariates)])])$sample.id

####LD pruning to get variant set
snp.dat <- data.frame(variant.id = seqGetData(gds, "variant.id"), chr = seqGetData(gds, "chromosome"), pos = seqGetData(gds, "position"))
snp.dat$chr_pos <- paste0(snp.dat$chr, ":", snp.dat$pos)
if(snpset.file=="null"){
  snpset <- snpgdsLDpruning(gds, sample.id=analysis.sample.id, method="corr", slide.max.bp=10e7, ld.threshold=sqrt(0.1))
  pruned <- unlist(snpset, use.names=FALSE)
  saveRDS(pruned, "pruned.rds")
  pruned.dat <- snp.dat[snp.dat$variant.id%in%pruned,]
  fwrite(pruned.dat[,c("chr", "pos")], "snpset.txt", quote=FALSE, row.names=FALSE, sep=',')
  seqSetFilter(gds, variant.id = pruned)
  seqExport(gds, pruned.gds.file)
}else{
  snpset.dat <- fread(snpset.file,stringsAsFactors=F,header=T,na.strings=c(NA,""))
  snpset.dat$chr_pos <- paste0(snpset.dat$chr, ":", snpset.dat$pos)
  snp.intersect.dat <- snp.dat[snp.dat$chr_pos%in% snpset.dat$chr_pos,]
  pruned <- unlist(snp.intersect.dat$variant.id)
  saveRDS(pruned, "pruned.rds")
  seqSetFilter(gds, variant.id = pruned)
  seqExport(gds, pruned.gds.file)
}
print(paste0("number of variants after pruning : ",length(pruned)))

date()
sink()
