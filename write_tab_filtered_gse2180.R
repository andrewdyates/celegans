load("~/c.elegans/GSE2180.MAR13.GOLDVAR.RData")
library("Biobase")

col.qq <- GSE2180.SCAN$genotype %in% c("N2", "ms")
max.expr <- apply(exprs(GSE2180.SCAN[,col.qq]), 1, max)
summary(max.expr)
gold.qq <- !is.na(featureData(GSE2180.SCAN)$gold)
summary(max.expr[gold.qq])
## remove probes with max expression less than 0.5 SCAN (lowest non-filtered gold is 0.645 max)

row.qq <- max.expr>=0.5
M <- exprs(GSE2180.SCAN[row.qq, col.qq])
dim(M) # 8995   61
write.table(M, file="~/c.elegans/GSE2180.SCAN.N2.ms.max_0.5.tab", sep="\t", row.names=T, col.names=T, quote=F)
save(M, file="~/c.elegans/GSE2180.SCAN.N2.ms.max_0.5.RData")