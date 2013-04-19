library(Biobase)
load("/nfs/01/osu6683/c.elegans/apr16.genelevel.depM.RData")
load("/nfs/01/osu6683/c.elegans/apr16.genelevel.exprs.RData")
# Local load of GSPLOM library
source("/nfs/01/osu6683/pymod/dependency_glyph_splom/lib.R")

# Gold
# ----------------------------------------
pdf("/nfs/01/osu6683/c.elegans/gold.gsplom.genes.apr17.pdf", width=15, height=15)
GOLD.R <- splom(D.expr.gold$CLS, D.expr.gold$DCOR, asGlyphs=T, useRaster=F, draw.labs=T, reorder=T, MAX=1, MIN=0.2)
dev.off()
GOLD.R$qq <- order.dendrogram(GOLD.R$Rhclust)
pdf("/nfs/01/osu6683/c.elegans/gold.dendro.genes.apr17.pdf", width=15, height=8)
plot(GOLD.R$Rhclust)
dev.off()


# Transcription factors
# ----------------------------------------
pdf("/nfs/01/osu6683/c.elegans/trans.gsplom.genes.apr17.pdf", width=40, height=40)
TRANS.R <- splom(D.expr.trans$CLS, D.expr.trans$DCOR, asGlyphs=T, useRaster=F, draw.labs=T, reorder=T, MAX=1, MIN=0.2)
dev.off()
pdf("/nfs/01/osu6683/c.elegans/trans.dendro.genes.apr17.pdf", width=50, height=8)
plot(TRANS.R$Rhclust)
dev.off()

# Plot As Compact
# ----------------------------------------
TRANS.R$qq <- order.dendrogram(TRANS.R$Rhclust)
pdf("/nfs/01/osu6683/c.elegans/trans.gsplom.genes.apr17.1px.pdf", width=40, height=40)
NULL.R <- splom(D.expr.trans$CLS[TRANS.R$qq,TRANS.R$qq], D.expr.trans$DCOR[TRANS.R$qq,TRANS.R$qq], asGlyphs=F, useRaster=F, draw.labs=T, reorder=F, MAX=1, MIN=0.2, grid=F)
dev.off()
save(TRANS.R, file="/nfs/01/osu6683/c.elegans/trans.genes.apr17.gsplom.RData")

# All expressed
# ----------------------------------------
png("/nfs/01/osu6683/c.elegans/all.gsplom.genes.apr17.png", width=5336, height=5336)
EXP.R <- splom(D.expr$CLS, D.expr$DCOR, asGlyphs=F, useRaster=T, draw.labs=F, reorder=T, MAX=1, MIN=0.3, grid=F)
dev.off()
save(EXP.R, file="/nfs/01/osu6683/c.elegans/all.genes.apr17.gsplom.RData")

pdf("/nfs/01/osu6683/c.elegans/all.dendro.genes.apr17.pdf", width=1000, height=30)
plot(EXP.R$Rhclust)
dev.off()