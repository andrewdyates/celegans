library("Biobase")
load("../celegans.apr8.expr.RData")

# NOTE: This will export phenotype attributes, too!
#write.table(E.gold, file="../celegans.apr8.gold.tab", sep="\t", col.names=NA, row.names=T)

write.table(exprs(E.gold), file="../celegans.apr8.gold.tab", sep="\t", col.names=NA, row.names=T, quote=F)
# assayData: 19 features, 61 samples 
write.table(exprs(E.expr.trans), file="../celegans.apr8.expr.trans.tab", sep="\t", col.names=NA, row.names=T, quote=F)
# assayData: 285 features, 61 samples 
write.table(exprs(E.expr), file="../celegans.apr8.expr.tab", sep="\t", col.names=NA, row.names=T, quote=F)
# assayData: 7063 features, 61 samples 
