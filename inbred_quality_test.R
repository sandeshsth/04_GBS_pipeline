#!/usr/bin/env Rscript

ARGV = commandArgs(trailingOnly=TRUE)

stats = read.table(file=ARGV[1], header=TRUE, check.names=FALSE)
sum(stats$inbredF>0.8)
score = round(stats$inbredF*100)
score[score<0]=0
sum(score>90)
qs = cbind(stats[,1:2], score)
colnames(qs) = c("CHROM", "POS", "QUALITYSCORE")
write.table(qs, file=paste0("inbred_quality_", ARGV[1]), quote=FALSE, sep="\t", row.names=FALSE)
