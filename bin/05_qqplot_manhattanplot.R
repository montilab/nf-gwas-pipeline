#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
summary.file <- args[1]
MAF <- as.numeric(args[2])
max_pval <- as.numeric(args[3])

sink('qqplot_manhattanplot.log', append=FALSE, split=TRUE)
date()
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(latex2exp))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(qqman))

####input data
dat <- fread(summary.file,header = T,stringsAsFactors = F)
dat <- dat %>%
  rename(pval = contains("pval"))
dat <- dat[which(!is.na(dat$pval) & dat$pval!=0),]
print(paste0("variants with valid results : ", dim(dat)[1]))

####qqplot
lambda_GC <- median(qchisq(1-dat$pval,1))/qchisq(0.5,1)
index <- 1:dim(dat)[1]/(dim(dat)[1]+1)

cat("\n####qq-plot starts\n")
png("qqplot.png")
plot(-log10(index), -log10(sort(dat$pval)), xlab=TeX(paste0("Expected ","$-log_{10}P$")), ylab=TeX(paste0("Observed ","$-log_{10}P$")), pch=16)
text(0.1*max(-log10(index)),0.9*(max(-log10(sort(dat$pval)))),
     TeX(paste0("$\\lambda = $", round(lambda_GC,4))))
abline(0,1,col="red")
dev.off()
cat("####qq-plot ends\n\n")

####manhattan plot
cat("\n####manhattan-plot starts\n")
png("manhattan_plot.png",width=1200,height=400)
manhattan(dat, chr="chr", bp="pos", snp="variant.id", p="pval" )
dev.off()
cat("####manhattan-plot ends\n\n")

####filter out MAF
dat <- dat[dat$caf > MAF & dat$caf < 1-MAF,]
print(paste0("variants with MAF > ", MAF, " : ", dim(dat)[1]))

####qqplot (filter out MAF)
lambda_GC <- median(qchisq(1-dat$pval,1))/qchisq(0.5,1)
index <- 1:dim(dat)[1]/(dim(dat)[1]+1)

cat("\n####qq-plot for MAF starts\n")
png(paste0("qqplot_MAF_",MAF,".png"))
plot(-log10(index), -log10(sort(dat$pval)), xlab=TeX(paste0("Expected ","$-log_{10}P$")), ylab=TeX(paste0("Observed ","$-log_{10}P$")), pch=16)
text(0.1*max(-log10(index)),0.9*(max(-log10(sort(dat$pval)))),
     TeX(paste0("$\\lambda = $", round(lambda_GC,4))))
abline(0,1,col="red")
dev.off()
cat("####qq-plot for MAF ends\n\n")

####manhattan plot (filter out MAF)
cat("\n####manhattan-plot for MAF starts\n")
png(paste0("manhattan_plot_MAF_",MAF,".png"),width=1200,height=400)
manhattan(dat, chr="chr", bp="pos", snp="variant.id", p="pval" )
dev.off()
cat("####manhattan-plot for MAF ends\n\n")

date()
sink()
