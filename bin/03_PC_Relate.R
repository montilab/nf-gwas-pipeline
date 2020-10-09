#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
gds.file <- args[1]
pheno.file <- args[2]
phenotypes <- args[3]
covariates <- unlist(strsplit(args[4], ","))
snpset.file <- args[5]
analysis.sample.id <- args[6]
annot <- args[7]
pruned <- args[8]
king <- args[9]
pcs <- args[10]
pc_df <- args[11]

sink("pc_relate.log", append=FALSE, split=TRUE)
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

####Read in pruned set
pruned <- readRDS(pruned)

####KING
king <- readRDS(king)
kingMat <- king$kinship
dimnames(kingMat) <- list(king$sample.id, king$sample.id)

####PC-AiR
pcs <- readRDS(pcs)

pc.df <- readRDS(pc_df)

####PC-Relate
seqSetFilter(seqData, variant.id=pruned)
iterator <- SeqVarBlockIterator(seqData, variantBlock=20000, verbose=FALSE)

cat("\n####pcrelate starts\n")
pcrel <- pcrelate(iterator, pcs=pcs$vectors[,1:2], sample.include=analysis.sample.id, training.set=pcs$unrels)
cat("####pcrelate ends\n\n")

seqResetFilter(seqData, verbose=FALSE)

kinship <- pcrel$kinBtwn

cat("\n####kinship plot starts\n")
png("kinship.png")
ggplot(kinship, aes(k0, kin)) +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color="grey") +
    geom_point(alpha=0.5) +
    ylab("kinship estimate") +
    theme_bw()
dev.off()
cat("####kinship plot ends\n\n")

####covariance matrix from pcrelate output
grm <- pcrelateToMatrix(pcrel, scaleKin=2)
saveRDS(grm,"grm.rds")

date()
sink()
