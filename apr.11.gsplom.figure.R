library("Biobase")
library("RColorBrewer")
library("gplots")
library("dynamicTreeCut")
load("../celegans.apr8.expr.RData")
load("../apr9.dcor.cls.RData")
heat.cols <- colorRampPalette(rev(brewer.pal(11,"RdBu")))
source("~/pymod/dependency_glyph_splom/lib.R")

# color gold standard probes
phaseI <- c("pal-1 174043_at", "pal-1 193341_at", "pal-1 193342_s_at")
phaseII.out <- c("tbx-9 190539_s_at", "scrt-1 192307_at")
phaseII <- c("cwn-1 188486_s_at", "elt-1 192655_s_at", "hnd-1 192707_at", "tbx-8 190477_at")
phaseIII <- c("elt-3 175801_at", "nob-1 175771_at", "nob-1 191468_s_at", "nhr-25 175663_at", "unc-120 188706_at", "hlh-1 193759_at")
phaseIV <- c("mab-21 188239_s_at")
phaseNA <- c("lin-26;lir-1 174037_at", "lin-26;lir-1 190309_at")
cols <- brewer.pal(6,"Set1")
classCols <- rep("#ffffff", dim(TRANS.CLS)[1])
pns <- rownames(TRANS.CLS)
classCols[pns %in% phaseI] <- cols[1] # red
classCols[pns %in% phaseII.out] <- cols[2] # blue 
classCols[pns %in% phaseII] <- cols[3] # green
classCols[pns %in% phaseIII] <- cols[4] # purple
classCols[pns %in% phaseIV] <- cols[5] # orange 
classCols[pns %in% phaseNA] <- cols[6] # yellow

## ----------------------------------------
## GSPLOM
## ----------------------------------------
pdf("trans.gsplom.apr11.complete.dcor1.pdf", width=60, height=60)
par(mar=c(10,5,5,5))
TRANS.R <- splom(TRANS.CLS, TRANS.DCOR, asGlyphs=T, draw.labs=T, reorder=T, grid=F, MAX=1, MIN=0.2, clust.meth="complete", main="All expressed transcription factors, rasterized, no grid, dcor weight 1", DCOR.weight=1, useRaster=T)
dev.off()
TRANS.R.qq <- order.dendrogram(TRANS.R$Rhclust)

w <- length(TRANS.R.qq)*10
png("trans.gsplom.apr11.complete.dcor1.png", width=w, height=w+5)
par(mar=c(10,5,5,5))
xR <- splom(TRANS.CLS[TRANS.R.qq,TRANS.R.qq], TRANS.DCOR[TRANS.R.qq,TRANS.R.qq], asGlyphs=F, draw.labs=T, reorder=F, grid=F, MAX=1, MIN=0.2, clust.meth="complete", main="All expressed transcription factors, rasterized, no grid, dcor weight 1", DCOR.weight=1, useRaster=T)
dev.off()

pdf("trans.gsplom.apr11.complete.dcor2.pdf", width=60, height=60)
par(mar=c(10,5,5,5))
TRANS.R2 <- splom(TRANS.CLS, TRANS.DCOR, asGlyphs=T, draw.labs=T, reorder=T, grid=F, MAX=1, MIN=0.2, clust.meth="complete", main="All expressed transcription factors, rasterized, no grid, dcor weight 2", DCOR.weight=2, useRaster=T)
dev.off()
TRANS.R2.qq <- order.dendrogram(TRANS.R2$Rhclust)


## ----------------------------------------
## PCC clustering (classCols)
## ----------------------------------------
D <- dist(TRANS.PCC)
H <- hclust(D, method="complete")
Z.d <- cutreeDynamic(H, minClusterSize=2, method="hybrid", distM=as.matrix(D), deepSplit=1)
Z.c <- cutree(H, h=7.5)
Zg.d <- cutreeDynamic(as.hclust(TRANS.R2$Rhclust), minClusterSize=2, method="hybrid", distM=as.matrix(TRANS.R2$D.row), deepSplit=1)
g <- 7.5/max(H$height)*max(as.hclust(TRANS.R2$Rhclust)$height) # 15.28
Zg.c.g <- cutree(as.hclust(TRANS.R2$Rhclust), h=g) # 93!!!
Zg.c <- cutree(as.hclust(TRANS.R2$Rhclust), k=15)

Z.d.cols <- colorRampPalette(brewer.pal(11,"Spectral"))(length(unique(Z.d)))[Z.d]
Z.c.cols <- colorRampPalette(brewer.pal(11,"Spectral"))(length(unique(Z.c)))[Z.c]
Zg.d.cols <- colorRampPalette(brewer.pal(11,"Spectral"))(length(unique(Zg.d)))[Zg.d]
Zg.c.cols <- colorRampPalette(brewer.pal(11,"Spectral"))(length(unique(Zg.c)))[Zg.c]


pdf("trans.pcc.heatmap.apr11.pdf", width=30, height=30)
# compare cut methods
R <- heatmap.2(TRANS.PCC, symbreaks=T, trace="none", col=heat.cols, Rowv=as.dendrogram(H), Colv=as.dendrogram(H), ColSideColors=Z.d.cols, RowSideColors=Z.c.cols, main="rows: 7.5 cut, cols: dynamic cut")
# compare height cut with baugh classes
R <- heatmap.2(TRANS.PCC, symbreaks=T, trace="none", col=heat.cols, Rowv=as.dendrogram(H), Colv=as.dendrogram(H), ColSideColors=classCols, RowSideColors=Z.c.cols, main="rows: 7.5 cut, cols: PHASE")
R <- heatmap.2(TRANS.PCC, symbreaks=T, trace="none", col=heat.cols, Rowv=as.dendrogram(H), Colv=as.dendrogram(H), ColSideColors=Zg.d.cols, RowSideColors=Z.d.cols, main="rows: dynamic cut, cols: GSPLOM dynamic cut")
R <- heatmap.2(TRANS.PCC, symbreaks=T, trace="none", col=heat.cols, Rowv=as.dendrogram(H), Colv=as.dendrogram(H), ColSideColors=Zg.c.cols, RowSideColors=Z.c.cols, main="rows: 7.5 cut, cols: GSPLOM k=15 cut")
R <- heatmap.2(TRANS.PCC, symbreaks=T, trace="none", col=heat.cols, Rowv=TRANS.R2$Rhclust, Colv=TRANS.R2$Rhclust, ColSideColors=Zg.c.cols, RowSideColors=Z.c.cols, main="GSPLOM dendro. rows: 7.5 cut, cols: GSPLOM k=15 cut")
dev.off()


pdf("gsplom.hclust.pdf", width=45, height=8)
par(mar=c(15,5,5,5))
plot(TRANS.R$Rhclust, main="dCor weight 1")
plot(TRANS.R2$Rhclust, main="dCor weight 2")
dev.off()

save(TRANS.R2, file="apr11.trans.gsplom.clust.RData")
## ----------------------------------------
## plot time series plots in clustering orders
## ----------------------------------------


