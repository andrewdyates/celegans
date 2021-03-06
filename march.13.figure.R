source("celegans.lib.R")
load("../GSE2180.MAR13.GOLDVAR.RData")
library(modeest)

## THESE SHOULD BE CLEANED UP SO AS NOT TO DEPEND ON gold.i

upc.maxs <- apply(exprs(GSE2180.UPC),1,max)
scan.stds <- apply(exprs(GSE2180.SCAN[upc.maxs>0.3,]),1,sd)
GSE2180.SCAN.b <- quantile(scan.stds,0.03)*2
geo.stds <- apply(exprs(GSE2180.GEO[upc.maxs>0.3,]),1,sd)
GSE2180.GEO.b <- quantile(geo.stds,0.03)*2

gold.i <- which(!is.na(featureData(GSE2180.SCAN)$gold))

## Plot differences between SCAN and GEO
## ----------------------------------------
plot.geo.scan.diff <- function(GEO, SCAN, name)
{
# scatterplots
pdf(paste0("gse2180.gold.geo.vs.scan.plots.",name,".pdf"))
for (i in gold.i) {
  x <- exprs(GEO)[i,]
  y <- exprs(SCAN)[i,]
  rho <- cor(x, y)
  d <- dcor(x, y)
  title <- paste(name, featureData(GEO)$gold[i], "PCC", format(rho,2), "dCOR", format(d,2))
  plot(x,y,main=title, xlab="GEO RMA norm", ylab="SCAN norm")
}
dev.off()
# time series
time <- SCAN$time
pdf(paste0("gse2180.time.mutant.gold.geo.vs.scan.plots.",name,".pdf"))
for (i in gold.i) {
  x <- exprs(GEO)[i,]
  y <- exprs(SCAN)[i,]
  rho <- cor(x, y)
  title <- paste(name, featureData(SCAN)$gold[i], "rho=", format(rho,digits=4))
  par(mar=c(5,4,4,5)+.1)
  plot(time, x, col="#cc0000", pch=as.numeric(SCAN$genotype), xlab="time", ylab="RMA Intensity", main=title)
  par(new=TRUE)
  plot(time, y, col="#000099", xaxt="n",yaxt="n",xlab="",ylab="",pch=as.numeric(SCAN$genotype))
  axis(4)
  mtext("SCAN Intensity",side=4,line=3)
  legend("topleft",col=c("#cc0000","#000099"),lty=1,legend=c("RMA","SCAN"))
  legend("topright",pch=1:4,legend=levels(SCAN$genotype))
}
dev.off()
# sploms
pdf(paste0("gse2180.geo.gold.pairs.",name,".pdf"), width=30, height=30)
pairs(t(exprs(SCAN[gold.i,])), labels=featureData(SCAN)$gold[gold.i], main=paste(name, "SCAN GSE2180 C.Elegans GSN"))
pairs(t(exprs(GEO[gold.i,])), labels=featureData(GEO)$gold[gold.i], main=paste(name, "GEO GSE2180 C.Elegans GSN"))
dev.off()
}

# Plot each mutant in time per gene
plot.time.series <- function(G,name)
{
time <- G$time
col.scale <- brewer.pal(4,"Set1")
cols <- col.scale[as.numeric(G$genotype)]
pdf(paste0("gse2180.time.mutant.gold.plots.",name,".pdf"))
for (i in gold.i) {
  y <- exprs(G)[i,]
  title <- paste(name, featureData(G)$gold[i])
  plot(time, y, col=cols, pch=20, main=title, ylab="Intensity")
  legend("topleft", pch=20, col=col.scale, legend=levels(G$genotype))
  plot(time, y, main=paste("No Labels", title), ylab="Intensity")
}
dev.off()
pdf(paste0("tiny.gse2180.time.mutant.gold.plots.",name,".pdf"),width=1,height=1)
for (i in gold.i) {
  y <- exprs(G)[i,]
  title <- paste(name, featureData(G)$gold[i])
  par(mar = c(0,0,0,0))
  plot(time, y, pch=20,ylab="", xlab="")
}
dev.off()
}

plot.compare.dep.matrix <- function(GEO.DCOR, GEO.PCC, SCAN.DCOR, SCAN.PCC, name) {
  pdf(paste0("dCOR.PCC.GEO.SCAN.hist.",name,".pdf"))
  hist(SCAN.DCOR-abs(SCAN.PCC), main=paste(name,"SCAN dCOR - SCAN |PCC|"))
  hist(SCAN.DCOR-GEO.DCOR, main=paste(name,"SCAN dCOR - GEO dCOR"))
  hist(GEO.DCOR-abs(GEO.PCC), main=paste(name,"GEO dCOR - GEO |PCC|"))
  hist(SCAN.DCOR-abs(GEO.PCC), main=paste(name,"SCAN dCOR - GEO |PCC|"))
  hist(abs(SCAN.PCC)-abs(GEO.PCC), main=paste(name,"SCAN |PCC| - GEO |PCC|"))
  dev.off()

  pdf(paste0("dCOR.PCC.SCAN.GEO.diff",name,".pdf"), width=20, height=20)
  heatmap.2(SCAN.DCOR - GEO.DCOR, main="SCAN DCOR - GEO DCOR", symm=T, col=heatmap_cols.near0, breaks=heatmap_breaks.near0, Colv=NULL, symm=T, trace='none')
  heatmap.2(abs(SCAN.PCC) - abs(GEO.PCC), main="SCAN |PCC| - GEO |PCC|", symm=T, col=heatmap_cols.near0, breaks=heatmap_breaks.near0, Colv=NULL, symm=T, trace='none')
  heatmap.2(SCAN.DCOR - abs(GEO.PCC), main="SCAN DCOR - GEO |PCC|", symm=T, col=heatmap_cols.near0, breaks=heatmap_breaks.near0, Colv=NULL, symm=T, trace='none')
  dev.off()
}

plot.dep.matrix <- function(R, name) {
  pdf(paste0("dependency.heatmaps.",name,".pdf"),width=10,height=10)
  heatmap.2(R$PCC[R$rowInd,R$rowInd], main=paste(name,"PCC"), symm=T, col=heatmap_cols, Rowv=NULL, dendrogram='none', trace='none', cexRow=0.6, cexCol=0.6)
  heatmap.2(R$DCOR[R$rowInd,R$rowInd], main=paste(name,"DCOR"), symm=T, col=heatmap_cols, breaks=heatmap_breaks, Rowv=NULL, dendrogram='none', trace='none', cexRow=0.6, cexCol=0.6)
  heatmap.2(R$COV[R$rowInd,R$rowInd], main=paste(name,"COV"), symm=T, col=heatmap_cols, Rowv=NULL, dendrogram='none', trace='none', cexRow=0.6, cexCol=0.6)
  heatmap.2(R$DCOV[R$rowInd,R$rowInd], main=paste(name,"dCOV"), symm=T, col=heatmap_cols, Rowv=NULL, dendrogram='none', trace='none', cexRow=0.6, cexCol=0.6)
  heatmap.2(R$DCOR[R$rowInd,R$rowInd] - abs(R$PCC[R$rowInd,R$rowInd]), main=paste(name,"dCOR - |PCC|"), symm=T, col=heatmap_cols, Rowv=NULL, dendrogram='none', trace='none', cexRow=0.6, cexCol=0.6)
  heatmap.2(R$DCOV[R$rowInd,R$rowInd] - abs(R$COV[R$rowInd,R$rowInd]), main=paste(name,"dCOV - |COV|"), symm=T, col=heatmap_cols, Rowv=NULL, dendrogram='none', trace='none', cexRow=0.6, cexCol=0.6)
  plot(R$Rhclust, main="Row Hcluster dCOR default clustering using dist+PCC")

  # Cluster on |PCC|, |dCOR|, add dendrograms.
  Z <- heatmap.2(R$DCOR, hclustfun=function(x) hclust(x,method="average"), symm=T, col=heatmap_cols, breaks=heatmap_breaks, trace='none', main="dCor average linkage euclidean distance")
  heatmap.2(R$PCC, hclustfun=function(x) hclust(x,method="average"), symm=T, col=heatmap_cols, trace='none', main="PCC average linkage euclidean distance")
  heatmap.2(R$SP, hclustfun=function(x) hclust(x,method="average"), symm=T, col=heatmap_cols, trace='none', main="SPEARMAN average linkage euclidean distance")
  heatmap.2(abs(R$PCC), hclustfun=function(x) hclust(x,method="average"), symm=T, col=heatmap_cols, breaks=heatmap_breaks, trace='none', main="abs(PCC) average linkage euclidean distance")
  heatmap.2(abs(R$SP), hclustfun=function(x) hclust(x,method="average"), symm=T, col=heatmap_cols, breaks=heatmap_breaks, trace='none', main="abs(SPEARMAN) average linkage euclidean distance")
  heatmap.2(matrix(0:9/9, nrow=5, ncol=2), col=heatmap_cols, breaks=heatmap_breaks, trace='none', main="force color scale")
  
  # dCor - |PCC|
  DIFF <- R$DCOR-abs(R$PCC)
  print(paste(name, range(DIFF)))
  x <- Z$rowInd
  heatmap.3(DIFF[x,x], reorder=F, MIN=-0.10, MAX=0.4, main="dCOR-|PCC| dCOR cluster order")
  heatmap.3(DIFF[x,x], reorder=F, MIN=-0.05, MAX=0.3, main="dCOR-|PCC| dCOR cluster order, color scale 2")
  heatmap.3(DIFF[x,x], reorder=F, MIN=min(DIFF), MAX=max(DIFF), main="dCOR-|PCC|, dCOR cluster order")
  DIFF.SP <- R$DCOR-abs(R$SP)
  print(paste(name, range(DIFF.SP), "SP"))
  heatmap.3(DIFF.SP[x,x], reorder=F, MIN=min(DIFF.SP), MAX=max(DIFF.SP), main="dCOR-|SP|, dCOR cluster order")
  dev.off()

  DIFF <- R$DCOR - abs(R$PCC)
  DIFF.SP <- R$DCOR - abs(R$SP)
  pdf(paste0("dependency.hists.",name,".pdf"))
  hist(R$PCC, main=paste(name, "PCC"))
  hist(R$SP, main=paste(name, "SPEARMAN"))
  hist(R$PCC[upper.tri(R$PCC)], main=paste(name, "PCC"))
  hist(R$DCOR, main=paste(name, "dCOR"))
  hist(R$DCOR[upper.tri(R$DCOR)], main=paste(name, "dCOR upper triangle"))
  hist(R$COV, main=paste(name, "COV"))
  hist(R$DCOV, main=paste(name, "dCOV"))
  hist(DIFF, main=paste(name, "dCOR - |PCC|"))
  hist(DIFF[upper.tri(DIFF)], main=paste(name, "dCOR - |PCC| upper triangle"))
  hist(DIFF.SP[upper.tri(DIFF.SP)], main=paste(name, "dCOR - |SP| upper triangle"))
  hist(R$DCOV - abs(R$COV), main=paste(name, "dCOV - |COV|"))
  dev.off()
}
  

label.dep.m <- function(D, labels) {
  rownames(D) <- labels
  colnames(D) <- labels
  D
}

calc.deps <- function(M, labels) {
  R <- list()
  R$DCOR <- label.dep.m(all.pairs.dcor(M), labels)
  R$DCOV <- label.dep.m(all.pairs.dcov(M), labels)
  R$PCC <- label.dep.m(cor(t(M)), labels)
  R$SP  <- label.dep.m(cor(t(M), method="spearman"), labels)
  R$COV <- label.dep.m(cov(t(M)), labels)

  # dcor hclust order
  D.DCOR.r <- as.dist(1-cor(t(R$DCOR), method="pearson")) + dist(R$DCOR)*2
  R$Rowv.dcor <- rowMeans(R$DCOR, na.rm = TRUE)
  R$Rhclust.dcor <- as.dendrogram(hclust(D.DCOR.r, method="average"))
  R$Rhclust.dcor <- reorder(R$Rhclust.dcor, R$Rowv.dcor)
  R$rowInd <- order.dendrogram(R$Rhclust.dcor)
  R
}

calc.plot.classes <- function(M, b, labels, name) {
  R <- list()
  rownames(M) <- labels
  pdf(paste0("stepfit.",name,".pdf"))
  R$steps <- all.steps(M)
  dev.off()
  pdf(paste0("boolclasses.",name,".pdf"))
  R$CLS <- all.pairs.cls(M, R$steps, b)
  dev.off()
  R
}

plot.pairs <- function(M, R, labels, name) {
  pdf(paste0("pairs.",name,".pdf"), width=20, height=20)
  pairs(t(M[R$rowInd,]), labels=labels[R$rowInd], main=name)
  dev.off()
}

plot.sploms <- function(CLS, D, name, MIN=0.18, scalar=10) {
  # TODO: return 
  D.DCOR.r <- as.dist(1-cor(t(D$DCOR), method="pearson")) + dist(D$DCOR)
  D.cls.r <- dist(CLS$CLS)
  Rowv <- rowMeans(D$DCOR, na.rm = TRUE)
  DCOR.weight <- 1
  Rhclust <- as.dendrogram(hclust(D.DCOR.r+D.cls.r, method="average"))
  Rhclust <- reorder(Rhclust, Rowv)
  CD.rowInd <- order.dendrogram(Rhclust)

  Rhclust.cls <- as.dendrogram(hclust(D.cls.r, method="average"))
  Rowv.cls <- apply(CLS$CLS, 1, function(x) unlist(mlv(as.vector(x), method="mfv")[1]))
  Rhclust.cls <- reorder(Rhclust.cls, Rowv.cls)
  C.rowInd <- order.dendrogram(Rhclust.cls)


  pdf(paste0("splom.",name,".pdf"), width=12, height=12)
  splom(CLS$CLS, D$DCOR, asGlyphs=T, MAX=1, MIN=0.18, reorder=F, main=paste0(name," Original Order, min dCOR=", MIN))
  splom(CLS$CLS[D$rowInd,D$rowInd], D$DCOR[D$rowInd,D$rowInd], asGlyphs=T, MAX=1, MIN=MIN, reorder=F, main=paste0(name, " dCOR Cluster Order, min dCOR=",MIN))
  splom(CLS$CLS[C.rowInd,C.rowInd], D$DCOR[C.rowInd,C.rowInd], asGlyphs=T, MAX=1, MIN=MIN, reorder=F, main=paste0(name, " CLS Cluster Order, min dCOR=",MIN))
  splom(CLS$CLS[CD.rowInd,CD.rowInd], D$DCOR[CD.rowInd,CD.rowInd], asGlyphs=T, MAX=1, MIN=MIN, reorder=F, main=paste0(name, " CLS+DCluster Cluster Order, min dCOR=",MIN))
  # GSPLOM with package clustering
  default.order.R <- splom(CLS$CLS, D$DCOR, asGlyphs=T, MAX=1, MIN=MIN, reorder=T, main=paste0(name, "Default order with CLS and DCOR, min dCOR=",MIN))
  # Different border options
  splom(CLS$CLS, D$DCOR, asGlyphs=T, MAX=1, MIN=MIN, reorder=T, main=paste0(name, "Default order with CLS lwd=3 and DCOR, min dCOR=",MIN), lwd=3)
  splom(CLS$CLS, D$DCOR, asGlyphs=T, MAX=1, MIN=MIN, reorder=T, main=paste0(name, "Default order with CLS lwd=3 and DCOR, min dCOR=",MIN), lwd=3, grid.col="#d3dfed")
  splom(CLS$CLS, D$DCOR, asGlyphs=T, MAX=1, MIN=MIN, reorder=T, main=paste0(name, "Default order with CLS lwd=5 and DCOR, min dCOR=",MIN), lwd=5)

  # no glyph form
  splom(CLS$CLS, D$DCOR, asGlyphs=F, MAX=1, MIN=MIN, reorder=T, main=paste0(name, "Default order with CLS and DCOR, min dCOR=",MIN))
  dev.off()

  pdf(paste0("splom.summary.",name,".pdf"), width=12, height=12)
  summary.plots(CLS$CLS,D$DCOR, sym=T)
  dev.off()

  width <- ncol(CLS$CLS)*scalar
  height <- nrow(CLS$CLS)*scalar
  png(paste0(name,'dft.raster.glyph.png'), width=width, height=height, units="px", bg="white")
  par(mar = rep(0, 4)) # set plot margins to 0 before drawing
  R <- splom(CLS$CLS, D$DCOR, asGlyphs=TRUE, useRaster=TRUE, draw.labs=F, reorder=T, MAX=1, MIN=MIN)
  dev.off()

  png(paste0(name,'dft.raster.pix.png'), width=width, height=height, units="px", bg="white")
  par(mar = rep(0, 4)) # set plot margins to 0 before drawing
  R <- splom(CLS$CLS, D$DCOR, asGlyphs=F, useRaster=TRUE, draw.labs=F, reorder=T, MAX=1, MIN=MIN)
  dev.off()

  pdf(paste0("dependency.dendros.",name,".pdf"),width=20,height=8)
  par(mar=c(10,0,5,0))
  plot(D$Rhclust.dcor, main=paste(name, "dCOR dendrogram"))
  plot(Rhclust, main=paste(name, "dCOR+CLS dendrogram"))
  plot(Rhclust.cls, main=paste(name, "CLS dendrogram"))
  plot(default.order.R$Rhclust, main=paste(name, "Default order row hclust"))
  plot(default.order.R$Chclust, main=paste(name, "Default order col hclust"))
  dev.off()

  RET <- list()
  RET$default.order.R <- default.order.R
  RET$CD.rowInd <- CD.rowInd
  RET
}

# --------------------
gold.labels <- featureData(GSE2180.GEO)$gold[gold.i] # same order as gold.i
analyze.and.plot <- function(GEO, SCAN, name) {
  R = list()
  plot.geo.scan.diff(GEO=GEO, SCAN=SCAN, name=name)
  plot.time.series(G=GEO, name=paste0("geo.",name))
  plot.time.series(G=SCAN, name=paste0("scan.",name))
  
  R$SCAN.D <- calc.deps(exprs(SCAN[gold.i,]), labels=gold.labels)
  R$GEO.D <- calc.deps(exprs(GEO[gold.i,]), labels=gold.labels)
  plot.dep.matrix(R$SCAN.D, name=paste0("SCAN.GOLD.",name))
  plot.dep.matrix(R$GEO.D, name=paste0("GEO.GOLD.",name))
  
  plot.pairs(M=exprs(SCAN[gold.i,]), R=R$SCAN.D, labels=gold.labels, name=paste0("SCAN.GOLD.",name))
  plot.pairs(M=exprs(GEO[gold.i,]), R=R$GEO.D, labels=gold.labels, name=paste0("GEO.GOLD.",name))
  
  R$SCAN.CLS <- calc.plot.classes(M=exprs(SCAN[gold.i,]), b=GSE2180.SCAN.b, labels=gold.labels, name=paste0("SCAN.GOLD.",name))
  R$GEO.CLS <- calc.plot.classes(M=exprs(GEO[gold.i,]), b=GSE2180.GEO.b, labels=gold.labels, name=paste0("GEO.GOLD.",name))
  
  R$SCAN.splom <- plot.sploms(CLS=R$SCAN.CLS, D=R$SCAN.D, name=paste0("SCAN.GOLD.",name))
  R$GEO.splom <- plot.sploms(CLS=R$GEO.CLS, D=R$GEO.D, name=paste0("GEO.GOLD.",name))

  # Scatterplot matrix of IEEE order
  o.row <- order.dendrogram(R$SCAN.splom$default.order.R$Rhclust)
  M <- exprs(SCAN[gold.i,])
  pdf(paste0("pairs.dtfgsplom.",name,".pdf"), width=20, height=20)
  pairs(t(M[o.row,]), main=paste(name,"Default GSPLOM order"))
  dev.off()
  R
}
## ==================================================

# 4mut
all.4mut <- analyze.and.plot(GEO=GSE2180.GEO, SCAN=GSE2180.SCAN, name="4mut")

# WT+ms
qq <- GSE2180.GEO$genotype %in% c('N2','ms')
wt.ms.2mut <- analyze.and.plot(GEO=GSE2180.GEO[,qq], SCAN=GSE2180.SCAN[,qq], name="WT+ms")

## ----------------------------------------
## HACK SCRIPT FOR PUBLICATION FIGURE GENERATION

R <- wt.ms.2mut$SCAN.D
DIFF <- R$DCOR - abs(R$PCC)
DIFF.SP <- R$DCOR - abs(R$SP)
name<-"wt.ms.2mut.scan"
pdf("wt.ms.2mut.scan.ieee.hists.pdf", width=6, height=6)
hist(R$DCOR[upper.tri(DIFF)], main=paste(name, "dCOR upper triangle"), breaks=seq(0,1,0.05))
hist(DIFF[upper.tri(DIFF)], main=paste(name, "dCOR - |PCC| upper triangle"), breaks=seq(-.2,.75,0.05))
hist(DIFF.SP[upper.tri(DIFF.SP)], main=paste(name, "dCOR - |SP| upper triangle"), breaks=seq(-.2,.75,0.05))
hist(R$DCOR[upper.tri(DIFF)], main=paste(name, "dCOR upper triangle"), breaks=seq(0,1,0.05), ylim=c(0,50))
hist(DIFF[upper.tri(DIFF)], main=paste(name, "dCOR - |PCC| upper triangle"), breaks=seq(-.2,.75,0.05), ylim=c(0,50))
hist(DIFF.SP[upper.tri(DIFF.SP)], main=paste(name, "dCOR - |SP| upper triangle"), breaks=seq(-.2,.75,0.05), ylim=c(0,50))

hist(R$PCC[upper.tri(R$PCC)], main=paste(name, "PCC upper triangle"), breaks=seq(-0.3,1,0.05), ylim=c(0,50))
hist(abs(R$PCC[upper.tri(R$PCC)]), main=paste(name, "abs(PCC) upper triangle"), breaks=seq(-0.3,1,0.05), ylim=c(0,50))
hist(R$SP[upper.tri(R$SP)], main=paste(name, "SP upper triangle"), breaks=seq(-0.3,1,0.05), ylim=c(0,50))
hist(abs(R$SP[upper.tri(R$SP)]), main=paste(name, "abs(SP) upper triangle"), breaks=seq(-0.3,1,0.05), ylim=c(0,50))
hist(R$DCOR[upper.tri(R$DCOR)], main=paste(name, "DCOR upper triangle"), breaks=seq(-0.3,1,0.05), ylim=c(0,50))

dev.off()

## dCOR, |PCC|, |SP| heatmaps in GSPLOM order
qq <- order.dendrogram(wt.ms.2mut$SCAN.splom$default.order.R$Rhclust)
R <- wt.ms.2mut$SCAN.D
DIFF.pcc <- R$DCOR - abs(R$PCC)
DIFF.sp  <- R$DCOR - abs(R$SP)
pdf("ieee_dependency_gsplom_order.pdf", width=10, height=10)
heatmap.3(R$DCOR[qq,qq], reorder=F, MIN=0, MAX=1, main="dCOR final GSPLOM order")
heatmap.3(DIFF.pcc[qq,qq], reorder=F, MIN=min(DIFF.pcc), MAX=max(DIFF.pcc), main="dCOR-|PCC| final GSPLOM order")
heatmap.3(DIFF.sp[qq,qq], reorder=F, MIN=min(DIFF.sp), MAX=max(DIFF.sp), main="dCOR-|SP| final GSPLOM order")
bdev.off()


## Bool glyphs only
pdf("ieee.cls.gsplom.no.dcorscaling.pdf", width=12, height=12)
# Cls distance only
par(mar=c(0,0,0,0))
Q <- splom.cls(wt.ms.2mut$SCAN.CLS$CLS, asGlyphs=T)
# GSPLOM order
par(mar=c(0,0,0,0))
Q2 <- splom.cls(wt.ms.2mut$SCAN.CLS$CLS[qq,qq], asGlyphs=T, reorder=F)
dev.off()
# WARNING: dendrograms do not include row/column names
Q.names <- rownames(wt.ms.2mut$SCAN.CLS$CLS)[order.dendrogram(Q$Rhclust)]


pdf("ieee.cls.distance.dendrograms.pdf", width=10, height=6)
plot(Q$Rhclust)
plot(Q$Chclust)
dev.off()

write.table(as.matrix(R$DCOR), file="gold.celegans.scan.dcor.tab", sep="\t")

pdf("ieee.splom.summary.pdf", width=8, height=6)
summary.plots(wt.ms.2mut$SCAN.CLS$CLS,wt.ms.2mut$SCAN.D$DCOR, sym=T)
dev.off()

## Export gold standard network as tab-delimited files.

qq <- order.dendrogram(wt.ms.2mut$SCAN.splom$default.order.R$Rhclust)
CLS <- wt.ms.2mut$SCAN.CLS$CLS[qq,qq]
DCOR <- wt.ms.2mut$SCAN.D$DCOR[qq,qq]
PCC <- wt.ms.2mut$SCAN.D$PCC[qq,qq]
qq.c <- GSE2180.GEO$genotype %in% c('N2','ms')

M <- exprs(GSE2180.SCAN[gold.i,qq.c])
M <- M[qq,] # in gold standard network order
time <- t(as.matrix(GSE2180.SCAN$time[qq.c]))
MM <- rbind(time, M)
rownames(MM)[1] <- "time(min)"
rownames(MM)[2:dim(MM)[1]] <- rownames(CLS)


# Export tables in excel friendly format
write.table(CLS, sep=",", file="gold.celegan.gse2180.cls.csv", col.names=NA, quote=F)
write.table(DCOR, sep=",", file="gold.celegan.gse2180.dcor.csv", col.names=NA, quote=F)
write.table(PCC, sep=",", file="gold.celegan.gse2180.pcc.csv", col.names=NA, quote=F)
write.table(MM, sep=",", file="gold.celegan.gse2180.M.csv", col.names=NA, quote=F)
gold.gsplom.dendro <- wt.ms.2mut$SCAN.splom$default.order.R$Rhclust
save(gold.gsplom.dendro, file="wt.ms.2mut.SCAN.gsplom.Rhclust.RData")
