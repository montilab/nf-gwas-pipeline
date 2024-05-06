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
data.plot <- c()
for(i in 1:22){
  print(i)
  temp <- dat[dat$chr == i & dat$pval < max_pval,]
  data.plot <- rbind(data.plot,temp[sort(temp$pos,index.return=T)[[2]],])
}

log10.pval <- -log10(data.plot$pval+0.) ## what goes on y-axis
color.length <- c(1,unlist(table(data.plot$chr)))
color.choice <- c(rep(c("darkred", "red"), 11),"darkred")

L.1 <- 0
L.2 <- 0
x.val <- c()
for(i in 1: (length(color.length)-1)){
  L.1 <- L.1 + color.length[i]
  L.2 <- L.2 + color.length[i+1]
  x.val <- c(x.val,(L.1+L.2)/2)
}

lab.chrom <- as.character(names(table(data.plot$chr)))
x<- c(1: length(log10.pval))
limit.bf <- min(log10.pval)

if(limit.bf==Inf){
	print("P-values not compatible")
	quit(status=0)
	}

L.1 <- 0
L.2 <- 0

png("manhattan_plot.png",width=1440,height=480,pointsize = 11)
par(mai=(c(0.65, 0.35, 0.35, 0.42) ))
plot(x,(log10.pval),ylim=c(limit.bf,max(log10.pval)),ylab="",axes=F,xlab="",cex=1.5,cex.lab=1.5,cex.axis=1.4)
axis(side=1, at=x.val, label=lab.chrom,cex.axis=1.4)
axis(side=2, cex.axis=1.5,line=-2.5)

for(i in 1: (length(color.length)-1)){
  L.1 <- L.1 + color.length[i]
  L.2 <- L.2 + color.length[i+1]
  points( x[L.1:L.2], log10.pval[ L.1 : L.2], col=color.choice[i]) 
}

abline(-log10(5e-8),0,col="grey",lty=1,cex=1.4)
abline(-log10(5e-6),0,col="grey",lty=2,cex=1.4)

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
data.plot <- c()
for(i in 1:22){
  print(i)
  temp <- dat[dat$chr == i & dat$pval < max_pval,]
  data.plot <- rbind(data.plot,temp[sort(temp$pos,index.return=T)[[2]],])
}

log10.pval <- -log10(data.plot$pval+0.) ## what goes on y-axis
color.length <- c(1,unlist(table(data.plot$chr)))
color.choice <- c(rep(c("darkred", "red"), 11),"darkred")

L.1 <- 0
L.2 <- 0
x.val <- c()
for(i in 1: (length(color.length)-1)){
  L.1 <- L.1 + color.length[i]
  L.2 <- L.2 + color.length[i+1]
  x.val <- c(x.val,(L.1+L.2)/2)
}

lab.chrom <- as.character(names(table(data.plot$chr)))
x<- c(1: length(log10.pval))
limit.bf <- min(log10.pval)
L.1 <- 0
L.2 <- 0

png(paste0("manhattan_plot_MAF_",MAF,".png"),width=1440,height=480,pointsize = 11)
par(mai=(c(0.65, 0.35, 0.35, 0.42) ))
plot(x,(log10.pval),ylim=c(limit.bf,max(log10.pval)),ylab="",axes=F,xlab="",cex=1.5,cex.lab=1.5,cex.axis=1.4)
axis(side=1, at=x.val, label=lab.chrom,cex.axis=1.4)
axis(side=2, cex.axis=1.5,line=-2.5)

for(i in 1: (length(color.length)-1)){
  L.1 <- L.1 + color.length[i]
  L.2 <- L.2 + color.length[i+1]
  points( x[L.1:L.2], log10.pval[ L.1 : L.2], col=color.choice[i]) 
}

abline(-log10(5e-8),0,col="grey",lty=1,cex=1.4)
abline(-log10(5e-6),0,col="grey",lty=2,cex=1.4)

dev.off()
cat("####manhattan-plot for MAF ends\n\n")

date()
sink()
