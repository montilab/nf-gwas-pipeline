#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
summary.file <- args[1]
qqplot <- "qqplot"

sink('qqplot.log', append=FALSE, split=TRUE)
date()
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(latex2exp))
suppressPackageStartupMessages(library(data.table))

####MAF
MAF = 0.01

####input data
dat <- fread(summary.file,header = T,stringsAsFactors = F)
dat <- as.tbl(dat)
dat <- dat %>%
  rename(pval = contains("pval"))
dat <- dat[which(!is.na(dat$pval)),]

colnames(dat)[colnames(dat)=="Est"] <- "beta"
colnames(dat)[colnames(dat)=="Score"] <- "beta"
colnames(dat)[colnames(dat)=="Est.SE"] <- "se"
colnames(dat)[colnames(dat)=="Score.SE"] <- "se"

####qqplot
lambda_GC <- median((dat$beta/dat$se)^2)/qchisq(0.5,1)
index <- 1:dim(dat)[1]/(dim(dat)[1]+1)

cat("\n####qq-plot starts\n")
png(paste(qqplot,".png",sep=""))
plot(-log10(index), -log10(sort(dat$pval)), xlab=TeX(paste0("Expected ","$-log_{10}P$")), ylab=TeX(paste0("Observed ","$-log_{10}P$")), pch=16)
text(0.1*max(-log10(index)),0.9*(max(-log10(sort(dat$pval)))),
     TeX(paste0("$\\lambda = $", round(lambda_GC,4))))
abline(0,1,col="red")
dev.off()
cat("####qq-plot ends\n\n")


date()
sink()
