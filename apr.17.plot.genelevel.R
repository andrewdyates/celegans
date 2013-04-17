library(Biobase)
load("../apr16.genelevel.depM.RData")
load("../apr16.genelevel.exprs.RData")
# Local load of GSPLOM library
source("~/pymod/dependency_glyph_splom/lib.R")

# GOLD
# ----------------------------------------
pdf("../gold.gsplom.genes.apr17.pdf", width=15, height=15)
GOLD.R <- splom(D.expr.gold$CLS, D.expr.gold$DCOR, asGlyphs=T, useRaster=F, draw.labs=T, reorder=T, MAX=1, MIN=0.2)
dev.off()
GOLD.R$qq <- order.dendrogram(TRANS.R$Rhclust)
pdf("../gold.dendro.genes.apr17.pdf", width=15, height=8)
plot(GOLD.R$Rhclust)
dev.off()

# Transcription factors
# ----------------------------------------
pdf("../trans.gsplom.genes.apr17.pdf", width=40, height=40)
TRANS.R <- splom(D.expr.trans$CLS, D.expr.trans$DCOR, asGlyphs=T, useRaster=F, draw.labs=T, reorder=T, MAX=1, MIN=0.2)
dev.off()
pdf("../trans.dendro.genes.apr17.pdf", width=50, height=8)
plot(TRANS.R$Rhclust)
dev.off()

# plot as compact
TRANS.R$qq <- order.dendrogram(TRANS.R$Rhclust)
pdf("../trans.gsplom.genes.apr17.1px.pdf", width=40, height=40)
NULL.R <- splom(D.expr.trans$CLS[TRANS.R$qq,TRANS.R$qq], D.expr.trans$DCOR[TRANS.R$qq,TRANS.R$qq], asGlyphs=F, useRaster=F, draw.labs=T, reorder=F, MAX=1, MIN=0.2, grid=F)
dev.off()
save(TRANS.R, file="../trans.genes.apr17.gsplom.RData")

# All expressed
# ----------------------------------------
png("../all.gsplom.genes.apr17.png", width=5336, height=5336)
EXP.R <- splom(D.expr$CLS, D.expr$DCOR, asGlyphs=F, useRaster=T, draw.labs=F, reorder=T, MAX=1, MIN=0.3, grid=F)
dev.off()
save(EXP.R, file="../all.genes.apr17.gsplom.RData")

pdf("../all.dendro.genes.apr17.pdf", width=1000, height=30)
plot(EXP.R$Rhclust)
dev.off()