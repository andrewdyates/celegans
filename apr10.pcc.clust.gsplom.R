load("../apr9.dcor.cls.RData")
library("RColorBrewer")
library("gplots")
library("dynamicTreeCut")
cols <- brewer.pal(11,"RdBu")

pdf("trans.pcc.heatmap.complete.pdf", width=30, height=30)
par(mar=c(10,10,0,0))
R <- heatmap.2(TRANS.PCC, symbreaks=T, trace="none", col=colorRampPalette(rev(cols)), hclustfun=function(x) hclust(x, method="complete"))
dev.off()

pdf("trans.pcc.heatmap.average.pdf", width=30, height=30)
par(mar=c(10,10,0,0))
R <- heatmap.2(TRANS.PCC, symbreaks=T, trace="none", col=colorRampPalette(rev(cols)), hclustfun=function(x) hclust(x, method="average"))
dev.off()

pdf("trans.pcc.heatmap.single.pdf", width=30, height=30)
par(mar=c(10,10,0,0))
R <- heatmap.2(TRANS.PCC, symbreaks=T, trace="none", col=colorRampPalette(rev(cols)), hclustfun=function(x) hclust(x, method="single"))
dev.off()

library("dynamicTreeCut")
D <- dist(TRANS.PCC)
H <- hclust(D, method="complete")
Z <- cutreeDynamic(H, minClusterSize=4, method="hybrid", distM=as.matrix(D), deepSplit=0)


cls.cols <- colorRampPalette(brewer.pal(11,"Spectral"))(length(unique(Z)))

pdf("trans.pcc.heatmap.complete.split.pdf", width=30, height=30)
R <- heatmap.2(TRANS.PCC, symbreaks=T, trace="none", col=colorRampPalette(rev(cols)), Rowv=as.dendrogram(H), Colv=as.dendrogram(H), ColSideColors=cls.cols[Z], RowSideColors=cls.cols[Z])
dev.off()

pdf("trans.pcc.heatmap.complete.height.plot.pdf")
plot(H$height)
dev.off()


Zcut <- cutree(H, h=7.8)
cls.cols <- colorRampPalette(brewer.pal(11,"Spectral"))(length(unique(Zcut))) # 13 clusters
pdf("trans.pcc.heatmap.complete.cut.pdf", width=30, height=30)
R <- heatmap.2(TRANS.PCC, symbreaks=T, trace="none", col=colorRampPalette(rev(cols)), Rowv=as.dendrogram(H), Colv=as.dendrogram(H), ColSideColors=cls.cols[Zcut], RowSideColors=cls.cols[Zcut], main="Cut at height 7.8. Complete linkage, euclidean distance.")
dev.off()


## For each of 13 clusters, compute GSPLOM
source("~/pymod/dependency_glyph_splom/lib.R")

# plot each gsplom, get order per gsplom
clust.idxs <- split(1:nrow(TRANS.PCC), Zcut)
Rs <- list()
i <- 1
for (clust in clust.idxs) {
  C <- TRANS.CLS[clust,clust]
  D <- TRANS.DCOR[clust,clust]
  title <- paste0("trans.pcc.",i,".gsplom.pdf"); i<-i+1
  pdf(title, width=20, height=20)
  R <- splom(C,D,asGlyphs=T,MIN=0.2,MAX=1)
  dev.off()
  Rs <- c(Rs,R)
}