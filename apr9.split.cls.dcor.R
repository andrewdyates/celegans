library("Biobase")
load("../mar20.gse2180.steps.cls.RData")
load("../celegans.apr8.expr.RData")
# load dCOR, PCC matrices computed from all expressed probes
load("../celegans.apr8.expr.pkl.DCOR.values.sig=None.abs=F.tab.RData")
DCOR <- as.matrix(M)
load("../celegans.apr8.expr.pkl.PEARSON.values.sig=None.abs=F.tab.RData")
PCC <- as.matrix(M)
# There is a mistake in dCor matrix; namely, it is square root of the actual value. Fix that.
DCOR <- DCOR ^ 2

# Plot all pairs dCor, PCC
qq <- sample.int(length(c(DCOR)), 100000)
png("DCOR_vs_PCC_all_exprs.png", width=2000, height=2000)
plot(c(DCOR)[qq], c(PCC)[qq])
dev.off()
# plot histograms
pdf("DCOR_PCC_hists.pdf", width=12, height=8)
hist(DCOR, breaks=100, main="DCOR, all expressed")
hist(PCC, breaks=200, main="PCC, all expressed")
hist(abs(PCC), breaks=100, main="|PCC|, all expressed")
hist(DCOR-abs(PCC), breaks=100, main="DCOR-|PCC|, all expressed")
dev.off()

## For gold, trans, exprs, select proper rows and columns
## ----------------------------------------
all(1:dim(PCC)[1] == match(rownames(PCC), rownames(DCOR))) # OK

# GOLD
qq.gold.D <- match(featureNames(featureData(E.gold)), rownames(PCC))
qq.gold.C <- match(featureNames(featureData(E.gold)), rownames(CLS))
# handle missing probe due to insufficient expression in gold set
exclude <- is.na(qq.gold.D)
qq.gold.D <- qq.gold.D[!exclude]
qq.gold.C <- qq.gold.C[!exclude]
GOLD.DCOR <- DCOR[qq.gold.D, qq.gold.D]
GOLD.CLS <- CLS[qq.gold.C, qq.gold.C]
rownames(GOLD.DCOR) <- featureData(E.gold[!exclude])$transfact.label
colnames(GOLD.DCOR) <- featureData(E.gold[!exclude])$transfact.label
rownames(GOLD.CLS) <- featureData(E.gold[!exclude])$transfact.label
colnames(GOLD.CLS) <- featureData(E.gold[!exclude])$transfact.label

# E.expr.trans (Expressed transcription factors)
qq.trans.D <- match(featureNames(featureData(E.expr.trans)), rownames(PCC))
qq.trans.C <- match(featureNames(featureData(E.expr.trans)), rownames(CLS))
TRANS.DCOR <- DCOR[qq.trans.D, qq.trans.D]
TRANS.CLS <- CLS[qq.trans.C, qq.trans.C]
rownames(TRANS.DCOR) <- featureData(E.expr.trans)$transfact.label
colnames(TRANS.DCOR) <- featureData(E.expr.trans)$transfact.label
rownames(TRANS.CLS) <- featureData(E.expr.trans)$transfact.label
colnames(TRANS.CLS) <- featureData(E.expr.trans)$transfact.label

# E.expr (All expressed probes)
qq.expr.D <- match(featureNames(featureData(E.expr)), rownames(PCC))
qq.expr.C <- match(featureNames(featureData(E.expr)), rownames(CLS))
EXPR.DCOR <- DCOR[qq.expr.D, qq.expr.D]
EXPR.CLS <- CLS[qq.expr.C, qq.expr.C]
rownames(EXPR.DCOR) <- featureData(E.expr)$transfact.label
colnames(EXPR.DCOR) <- featureData(E.expr)$transfact.label
rownames(EXPR.CLS) <- featureData(E.expr)$transfact.label
colnames(EXPR.CLS) <- featureData(E.expr)$transfact.label


# saved workspace. open to restore...
save(GOLD.DCOR, GOLD.CLS, TRANS.DCOR, TRANS.CLS, EXPR.DCOR, EXPR.CLS, file="../apr9.dcor.cls.RData")