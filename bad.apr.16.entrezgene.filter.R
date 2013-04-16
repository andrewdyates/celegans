# THIS IS BAD. DON'T USE IT
# ------------------------------
library(Biobase)
library(celegansceentrezg.db)
library(celegansceentrezgprobe)
load("../GSE2180.ALL.RData")
load("../celegans.apr8.expr.RData")

AttrT <- read.table("../GSE2180_GPL200.samples.tab", sep="\t", header=TRUE, row.names=1, as.is=T)
AttrT <- as.data.frame(t(AttrT))
AttrT$time <- as.numeric(as.character(AttrT$"n:Sample_time"))
AttrT$genotype <- AttrT$"n:Sample_genotype"

colnames(exprs(GSE2180.ALL.SCAN)) <- sub(".CEL.gz","",colnames(exprs(GSE2180.ALL.SCAN)))
colnames(exprs(GSE2180.ALL.UPC)) <- sub(".CEL.gz","",colnames(exprs(GSE2180.ALL.UPC)))

qq.c <- match(colnames(exprs(GSE2180.ALL.SCAN)), rownames(AttrT))
pData(phenoData(GSE2180.ALL.SCAN)) <- AttrT[qq.c,]

ids <- rownames(exprs(GSE2180.ALL.SCAN))
genes <- mget(ids, celegansceentrezgSYMBOL, ifnotfound = NA)
featureData(GSE2180.ALL.SCAN)$gene <- genes

gold.genes <- c('pal-1', 'tbx-8', 'tbx-9', 'elt-1', 'hnd-1', 'scrt-1', 'cwn-1', 'unc-120', 'hlh-1', 'nob-1', 'elt-3', 'nhr-25', 'mab-21', 'lin-26', 'vab-7')