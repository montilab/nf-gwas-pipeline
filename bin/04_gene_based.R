#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
gds.file <- args[1]
annot <- args[2]
nullmod <- args[3]
max_maf <- args[4]
method <- args[5]
result.file1 <- args[6]
result.file2 <- args[7]
log.file <- args[8]

sink(log.file, append=FALSE, split=TRUE)
date()
suppressPackageStartupMessages(library(SeqArray))
suppressPackageStartupMessages(library(GENESIS))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(SeqVarTools))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))

####Open GDS
gds <- seqOpen(gds.file)

####Create a SeqVarData object
annot <- readRDS(annot)
seqData <- SeqVarData(gds, sampleData=annot)

####Null model
nullmod <- readRDS(nullmod)

####Aggregate test
gr <- granges(gds)

####find variants that overlap with each gene
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
gr <- renameSeqlevels(gr, paste0("chr", seqlevels(gr)))
ts <- transcriptsByOverlaps(txdb, gr, columns=c("GENEID"))
genes.gr <- genes(txdb)

o <- findOverlaps(ts, genes.gr)

if(length(o)>0){
	gr2 <- split(genes.gr[subjectHits(o)], 1:length(o))

	####define genes
	genes <- unique(unlist(gr2))
	genes <- renameSeqlevels(genes, sub("chr", "", seqlevels(genes)))

	####create an iterator where each successive unit is a different gene
	iterator <- SeqVarRangeIterator(seqData, variantRanges = genes, verbose=FALSE)

	####do a burden test on the rare variants in each gene
	cat("\n####assocTestAggregate starts\n")
	assoc <- assocTestAggregate(iterator, nullmod, AF.max=max_maf, 
								test="Burden",
                                verbose=T)
	cat("####assocTestAggregate ends\n\n")

	out <- assoc$results
	out <- cbind(as.data.frame(genes), out)
	colnames(out)[1] <- "chr" 

	write.csv(out, file = result.file1, row.names=FALSE)
	saveRDS(assoc, result.file2)

}

if(length(o)==0){
	print("No genes identified")	
}

date()
sink()
