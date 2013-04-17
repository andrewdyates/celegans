library("Biobase")
load("../mar20.gse2180.steps.cls.RData")
load("../apr16.genelevel.exprs.RData")
# load dCOR, PCC matrices computed from all expressed probes
load("../celegans.apr8.expr.pkl.DCOR.values.sig=None.abs=F.tab.RData")
DCOR <- as.matrix(M)
load("../celegans.apr8.expr.pkl.PEARSON.values.sig=None.abs=F.tab.RData")
PCC <- as.matrix(M)
# There is a mistake in dCor matrix; namely, it is square root of the actual value. Fix that.
DCOR <- DCOR ^ 2


D.expr <- list()
qq <- match(rownames(exprs(Eg.expr)), rownames(CLS))
D.expr$CLS <- CLS[qq,qq]
qq <- match(rownames(exprs(Eg.expr)), rownames(DCOR))
D.expr$DCOR <- DCOR[qq,qq]
qq <- match(rownames(exprs(Eg.expr)), rownames(PCC))
D.expr$PCC <- PCC[qq,qq]

D.expr.gold <- list()
qq <- match(rownames(exprs(Eg.expr.gold)), rownames(CLS))
D.expr.gold$CLS <- CLS[qq,qq]
qq <- match(rownames(exprs(Eg.expr.gold)), rownames(DCOR))
D.expr.gold$DCOR <- DCOR[qq,qq]
qq <- match(rownames(exprs(Eg.expr.gold)), rownames(PCC))
D.expr.gold$PCC <- PCC[qq,qq]

D.expr.trans <- list()
qq <- match(rownames(exprs(Eg.expr.trans)), rownames(CLS))
D.expr.trans$CLS <- CLS[qq,qq]
qq <- match(rownames(exprs(Eg.expr.trans)), rownames(DCOR))
D.expr.trans$DCOR <- DCOR[qq,qq]
qq <- match(rownames(exprs(Eg.expr.trans)), rownames(PCC))
D.expr.trans$PCC <- PCC[qq,qq]

save(D.expr, D.expr.gold, D.expr.trans, file="../apr16.genelevel.depM.RData")