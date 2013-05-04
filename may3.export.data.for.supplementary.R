load("../apr16.genelevel.exprs.RData")
# we want a nice table of gene expression and gene names

# Gold network genes only
# ----------------------------------------
write.table(Eg.expr.gold, file="../raw.may3.Eg.expr.gold.celegans.csv", sep=",", col.names=NA, row.names=T)
Z <- pData(featureData(Eg.expr.gold))
write.table(as.matrix(Z), file="../raw.may3.Eg.expr.gold.celegans.featureData.tab", sep="\t", col.names=NA, row.names=T)
# for regular people
M <- exprs(Eg.expr.gold)
rownames(M) <- featureData(Eg.expr.gold)$sym
write.table(M, file="../nice.may3.Eg.expr.gold.celegans.csv", sep=",", col.names=NA, row.names=T)

# Same thing for all transcription factors
# ----------------------------------------
write.table(Eg.expr.trans, file="../raw.may3.Eg.expr.alltrans.celegans.csv", sep=",", col.names=NA, row.names=T)
Z <- pData(featureData(Eg.expr.trans))
write.table(as.matrix(Z), file="../raw.may3.Eg.expr.alltrans.celegans.featureData.tab", sep="\t", col.names=NA, row.names=T)
M <- exprs(Eg.expr.trans)
rownames(M) <- featureData(Eg.expr.trans)$sym
write.table(M, file="../nice.may3.Eg.expr.alltrans.celegans.csv", sep=",", col.names=NA, row.names=T)


#Eg.expr.trans
