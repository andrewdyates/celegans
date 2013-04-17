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


get.dep.M <- function(Eg, CLS, DCOR, PCC) {
  R <- list()
  R$rownames <- rownames(exprs(Eg))
  R$syms <- featureData(Eg)$sym
  qq <- match(R$rownames, rownames(CLS))
  R$CLS <- CLS[qq,qq]
  rownames(R$CLS) <- R$syms
  colnames(R$CLS) <- R$syms
  qq <- match(R$rownames, rownames(DCOR))
  R$DCOR <- DCOR[qq,qq]
  rownames(R$DCOR) <- R$syms
  colnames(R$DCOR) <- R$syms
  qq <- match(R$rownames, rownames(PCC))
  R$PCC <- PCC[qq,qq]
  rownames(R$PCC) <- R$syms
  colnames(R$PCC) <- R$syms
  R
}

D.expr <- get.dep.M(Eg.expr, CLS, DCOR, PCC)
D.expr.trans <- get.dep.M(Eg.expr.trans, CLS, DCOR, PCC)
D.expr.gold <- get.dep.M(Eg.expr.gold, CLS, DCOR, PCC)

save(D.expr, D.expr.gold, D.expr.trans, file="../apr16.genelevel.depM.RData")