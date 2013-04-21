library(Biobase)
load("../trans.genes.apr17.gsplom.RData")
load("../apr16.genelevel.depM.RData")
load("../apr16.genelevel.exprs.RData")
source("~/Dropbox/biostat/git_repos/dependency_glyph_splom/lib.R")

# TRANS
# GSPLOM with complete linkage, no glyphs on TRANS
w <- dim(D.expr.trans$CLS)[1]
pdf("~/Desktop/expr.trans.gsplom.pdf", width=50, height=50)
R.TRANS <- splom(D.expr.trans$CLS, D.expr.trans$DCOR, asGlyphs=F, pad=F, grid=F, useRaster=F, clust.meth="complete", MIN=0.25, MAX=1, reorder=T)
dev.off()
save(R.TRANS, file="../apr.19.gsplom.R.RData")
pdf("~/Desktop/expr.trans.gsplom.dendro.pdf", width=60, height=8)
plot(R.TRANS$Rhclust)
dev.off()

# export GOLD dCor and PCC as .tab
write.table(D.expr.gold$CLS, file="../D.expr.gold.CLS.apr.19.tab", quote=F, sep="\t")
write.table(D.expr.gold$DCOR, file="../D.expr.gold.DCOR.apr.19.tab", quote=F, sep="\t")

# implement weak coupling


# run weak coupling on GOLD, export as .tab

