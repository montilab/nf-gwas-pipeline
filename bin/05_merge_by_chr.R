#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
result.file <- args[1]
caf.file <- args[2] 
model <- args[3]
combine.file <- args[4]
log.file <- args[5]

sink(log.file, append=FALSE, split=TRUE)
date()
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))

result.dat <- fread(result.file, stringsAsFactors = F, header = T)
caf.dat <- fread(caf.file, stringsAsFactors = F, header = T)

combine.dat <- left_join(result.dat, caf.dat, by=c("chr","pos"))
print(paste0("total variants analyzed : ", dim(combine.dat)[1]))
maf <- ifelse(combine.dat$caf>0.5, 1-combine.dat$caf, combine.dat$caf)
out.dat <- combine.dat[maf*2*combine.dat$N>0.99,]
print(paste0("variants with at least 1 copy : ", dim(out.dat)[1]))

if (model == "logistic") {
	case.maf <- ifelse(out.dat$case.caf>0.5, 1-out.dat$case.caf, out.dat$case.caf)
	control.maf <- ifelse(out.dat$control.caf>0.5, 1-out.dat$control.caf, out.dat$control.caf)

	out.dat <- out.dat[(case.maf*2*out.dat$n.case>0.99)|(control.maf*2*out.dat$n.control>0.99),]
	print(paste0("variants with at least 1 copy in either cases or controls : ", dim(out.dat)[1]))

}

fwrite(out.dat, combine.file, row.names=F)

date()
sink()
