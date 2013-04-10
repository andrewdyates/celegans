library("Biobase")
load("../celegans.apr8.expr.RData")
#GOLD.DCOR, GOLD.CLS, TRANS.DCOR, TRANS.CLS, EXPR.DCOR, EXPR.CLS
load("../apr9.dcor.cls.RData")
# this should be rerouted to a library import
source("~/pymod/dependency_glyph_splom/lib.R")

# gold
pdf("gold.gsplom.apr9.pdf", width=15, height=15)
GOLD.R <- splom(GOLD.CLS, GOLD.DCOR, asGlyphs=T, useRaster=F, draw.labs=T, reorder=T, MAX=1, MIN=0.2)
dev.off()

pdf("trans.gsplom.apr9.avg.pdf", width=60, height=60)
par(mar=c(10,5,5,5))
TRANS.R.A <- splom(TRANS.CLS, TRANS.DCOR, asGlyphs=T, useRaster=F, draw.labs=T, reorder=T, grid=F, MAX=1, MIN=0.2, clust.meth="average", main="All transcription factors, avg linkage")
dev.off()

pdf("trans.gsplom.apr9.complete.pdf", width=60, height=60)
par(mar=c(10,5,5,5))
TRANS.R.C <- splom(TRANS.CLS, TRANS.DCOR, asGlyphs=T, useRaster=F, draw.labs=T, reorder=T, grid=F, MAX=1, MIN=0.2, clust.meth="complete", main="All transcription factors, complete linkage")
dev.off()

pdf("trans.gsplom.apr9.single.pdf", width=60, height=60)
par(mar=c(10,5,5,5))
TRANS.R.S <- splom(TRANS.CLS, TRANS.DCOR, asGlyphs=T, useRaster=F, draw.labs=T, reorder=T, grid=F, MAX=1, MIN=0.2, clust.meth="single", main="All transcription factors, single linkage")
dev.off()

# de-weight dCor
pdf("trans.gsplom.apr9.avg.dcor1.pdf", width=60, height=60)
par(mar=c(10,5,5,5))
TRANS.R.A.1 <- splom(TRANS.CLS, TRANS.DCOR, asGlyphs=T, useRaster=F, draw.labs=T, reorder=T, grid=F, MAX=1, MIN=0.2, clust.meth="average", main="All transcription factors, avg linkage", DCOR.weight=1)
dev.off()

pdf("trans.gsplom.apr9.complete.dcor1.pdf", width=60, height=60)
par(mar=c(10,5,5,5))
TRANS.R.C.1 <- splom(TRANS.CLS, TRANS.DCOR, asGlyphs=T, useRaster=F, draw.labs=T, reorder=T, grid=F, MAX=1, MIN=0.2, clust.meth="complete", main="All transcription factors, complete linkage", DCOR.weight=1)
dev.off()

pdf("trans.gsplom.apr9.single.dcor1.pdf", width=60, height=60)
par(mar=c(10,5,5,5))
TRANS.R.S.1 <- splom(TRANS.CLS, TRANS.DCOR, asGlyphs=T, useRaster=F, draw.labs=T, reorder=T, grid=F, MAX=1, MIN=0.2, clust.meth="single", main="All transcription factors, single linkage", DCOR.weight=1)
dev.off()



TRANS.gsplom.qq <- order.dendrogram(TRANS.R$Rhclust)

pdf("trans.gsplom.apr9.pdf", width=60, height=60)
R <- splom(TRANS.CLS[TRANS.gsplom.qq,TRANS.gsplom.qq], TRANS.DCOR[TRANS.gsplom.qq,TRANS.gsplom.qq], asGlyphs=T, useRaster=F, draw.labs=T, reorder=F, MAX=1, MIN=0.2)
dev.off()
# no border
pdf("trans.gsplom.apr9.noborder.pdf", width=60, height=60)
R <- splom(TRANS.CLS[TRANS.gsplom.qq,TRANS.gsplom.qq], TRANS.DCOR[TRANS.gsplom.qq,TRANS.gsplom.qq], asGlyphs=T, useRaster=F, draw.labs=T, reorder=F, MAX=1, MIN=0.2, grid=F)
dev.off()
# raster (faster, smaller)
pdf("trans.gsplom.apr9.noborder.raster.pdf", width=60, height=60)
R <- splom(TRANS.CLS[TRANS.gsplom.qq,TRANS.gsplom.qq], TRANS.DCOR[TRANS.gsplom.qq,TRANS.gsplom.qq], asGlyphs=T, useRaster=T, draw.labs=T, reorder=F, MAX=1, MIN=0.2, grid=F)
dev.off()

w <- dim(TRANS.CLS)[1]
png("trans.gsplom.apr9.png", width=w, height=w)
par(mar=c(0,0,0,0))
R <- splom(TRANS.CLS[TRANS.gsplom.qq,TRANS.gsplom.qq], TRANS.DCOR[TRANS.gsplom.qq,TRANS.gsplom.qq], asGlyphs=F, useRaster=F, draw.labs=F, reorder=F, MAX=1, MIN=0.2, lwd=0)
dev.off()

# dCor heatmaps
pdf("trans.dcor.heatmap.apr9.pdf", width=60, height=60)
R <- heatmap.3(TRANS.DCOR[TRANS.gsplom.qq,TRANS.gsplom.qq], MIN=0.2, MAX=1, symm=T, reorder=F, main="GSPLOM order")
TRANS.R.HEAT <- heatmap.3(TRANS.DCOR, MIN=0.2, MAX=1, symm=T, reorder=T, main="dCor order")
dev.off()


pdf("trans.dcor.heatmap.apr9.pdf", width=60, height=60)
R <- heatmap.3(TRANS.DCOR[TRANS.gsplom.qq,TRANS.gsplom.qq], MIN=0.2, MAX=1, symm=T, reorder=F, main="GSPLOM order")
TRANS.R.HEAT <- heatmap.3(TRANS.DCOR, MIN=0.2, MAX=1, symm=T, reorder=T, main="dCor order")
dev.off()



pdf("trans.gsplom.summary.apr9.pdf", width=10, height=10)
summary.plots(TRANS.CLS, TRANS.DCOR, sym=T)
dev.off()

pdf("trans.gsplom.dendro.pdf", width=50, height=15)
par(mar=c(10,5,5,5))
plot(TRANS.R$Rhclust)
dev.off()

save(TRANS.R, file="apr9.TRANS.GSPLOM.results.RData")

TRANS.max <- apply(upper.tri(TRANS.DCOR),1,max)
TRANS.medians <- apply(TRANS.DCOR,1,median)
TRANS.means <- apply(TRANS.DCOR,1,mean)

