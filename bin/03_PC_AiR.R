å#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
gds.file <- args[1]
pheno.file <- args[2]
phenotypes <- args[3]
covariates <- unlist(strsplit(args[4], ","))
snpset.file <- args[5]

sink("pc_air.log", append=FALSE, split=TRUE)
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
print(paste0("number of individuals to analyze : ",length(analysis.sample.id)))
saveRDS(analysis.sample.id, "analysis.sample.id.rds")

metadata <- data.frame(labelDescription=colnames(annot),row.names=names(annot))
annot <- AnnotatedDataFrame(annot, metadata)
saveRDS(annot, "annot.rds")

all.equal(annot$sample.id,seqGetData(gds, "sample.id"))
seqData <- SeqVarData(gds, sampleData=annot)

####LD pruning to get variant set
snp.dat <- data.frame(variant.id = seqGetData(gds, "variant.id"), chr = seqGetData(gds, "chromosome"), pos = seqGetData(gds, "position"))
snp.dat$chr_pos <- paste0(snp.dat$chr, ":", snp.dat$pos)
if(snpset.file=="null"){
  snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=10e7, ld.threshold=sqrt(0.1))
  pruned <- unlist(snpset, use.names=FALSE)
  saveRDS(pruned, "pruned.rds")
  pruned.dat <- snp.dat[snp.dat$variant.id%in%pruned,]
  fwrite(pruned.dat[,c("chr", "pos")], "snpset.txt", quote=FALSE, row.names=FALSE, sep=',')
}else{
  snpset.dat <- fread(snpset.file,stringsAsFactors=F,header=T,na.strings=c(NA,""))
  snpset.dat$chr_pos <- paste0(snpset.dat$chr, ":", snpset.dat$pos)
  snp.intersect.dat <- snp.dat[snp.dat$chr_pos%in% snpset.dat$chr_pos,]
  pruned <- unlist(snp.intersect.dat$variant.id)
  saveRDS(pruned, "pruned.rds")
}
print(paste0("number of variants after pruning : ",length(pruned)))

####KING
cat("\n####snpgdsIBDKING starts\n")
king <- snpgdsIBDKING(gds, sample.id=analysis.sample.id, snp.id=pruned, verbose=T)
kingMat <- king$kinship
dimnames(kingMat) <- list(king$sample.id, king$sample.id)
saveRDS(king,"king.rds")
cat("####snpgdsIBDKING ends\n\n")

####PC-AiR
cat("\n####pcair starts\n")
pcs <- pcair(seqData, kinobj=kingMat, kin.thresh=2^(-7/2),
                      divobj=kingMat, div.thresh=-2^(-7/2),
             sample.include=analysis.sample.id,
             snp.include=pruned,
             verbose=T)
saveRDS(pcs,"pcs.rds")
cat("####pcair ends\n\n")

pc.df <- as.data.frame(pcs$vectors)
names(pc.df) <- paste0("PC", 1:ncol(pcs$vectors))
pc.df$sample.id <- row.names(pcs$vectors)
pc.df <- left_join(pData(annot), pc.df, by="sample.id")
saveRDS(pc.df,"pc.df.rds")

cat("\n####PCA plot starts\n")
png("PC1vsPC2.png")
ggplot(pc.df, aes(PC1, PC2)) + geom_point()
dev.off()

png("PC3vsPC4.png")
ggplot(pc.df, aes(PC3, PC4)) + geom_point()
dev.off()
cat("####PCA pot ends\n\n")

date()
sink()
