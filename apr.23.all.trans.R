library(Biobase)
load("../trans.genes.apr17.gsplom.RData")
load("../apr16.genelevel.depM.RData")
load("../apr16.genelevel.exprs.RData")
source("~/Dropbox/biostat/git_repos/dependency_glyph_splom/lib.R")
source("all.pairs.weak.R")

D.expr.trans$WEAK <- all.pairs.weak(exprs(Eg.expr.trans))
rownames(D.expr.trans$WEAK) <- featureData(Eg.expr.trans)$sym
colnames(D.expr.trans$WEAK) <- featureData(Eg.expr.trans)$sym

WEAK <- D.expr.trans$WEAK
save(WEAK, file="../trans.weak.RData")

all(rownames(D.expr.trans$CLS) == rownames(D.expr.trans$DCOR))
all(rownames(D.expr.trans$CLS) == rownames(WEAK))
all(rownames(D.expr.trans$CLS) == colnames(WEAK))

write.table(WEAK, file="../trans.weak.tab", sep="\t", quote=F, col.names=NA, row.names=T)
write.table(D.expr.trans$CLS, file="../trans.cls.tab", sep="\t", quote=F, col.names=NA, row.names=T)
write.table(D.expr.trans$DCOR, file="../trans.dcor.tab", sep="\t", quote=F, col.names=NA, row.names=T)
