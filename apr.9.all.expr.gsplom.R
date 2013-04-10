library("Biobase")
load("../celegans.apr8.expr.RData")
#GOLD.DCOR, GOLD.CLS, TRANS.DCOR, TRANS.CLS, EXPR.DCOR, EXPR.CLS
load("../apr9.dcor.cls.RData")
# this should be rerouted to a library import
source("~/pymod/dependency_glyph_splom/lib.R")


pdf("gold.gsplom.apr9.pdf", width=60, height=60)
par(mar=c(0,0,0,0))
R <- splom(EXPR.CLS, EXPR.DCOR, asGlyphs=T, useRaster=T, draw.labs=F, reorder=T, MAX=1, MIN=0.2)
dev.off()

png(paste0(name,'dft.raster.glyph.png'), width=width, height=height, units="px", bg="white")
par(mar = rep(0, 4)) # set plot margins to 0 before drawing
EXPR.R <- splom(EXPR.CLS, EXPR.DCOR, asGlyphs=TRUE, useRaster=T, draw.labs=F, reorder=T, MAX=1, MIN=0.2)
dev.off()