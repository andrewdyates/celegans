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
}
dev.off()
}

plot.dep.matrix <- function(D, name) {
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
  heatmap.2(SCAN.DCOR - GEO.DCOR), main="SCAN DCOR - GEO DCOR", symm=T, col=heatmap_cols.near0, breaks=heatmap_breaks.near0, Colv=NULL, symm=T, trace='none')
  heatmap.2(abs(SCAN.PCC) - abs(GEO.PCC)), main="SCAN |PCC| - GEO |PCC|", symm=T, col=heatmap_cols.near0, breaks=heatmap_breaks.near0, Colv=NULL, symm=T, trace='none')
  heatmap.2(SCAN.DCOR - abs(GEO.PCC)), main="SCAN DCOR - GEO |PCC|", symm=T, col=heatmap_cols.near0, breaks=heatmap_breaks.near0, Colv=NULL, symm=T, trace='none')
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
  dev.off()
  pdf(paste0("dependency.hists.",name,".pdf"))
  hist(R$PCC, main=paste(name, "PCC"))
  hist(R$DCOR, main=paste(name, "dCOR"))
  hist(R$COV, main=paste(name, "COV"))
  hist(R$DCOV, main=paste(name, "dCOV"))
  hist(R$DCOR - abs(R$PCC), main=paste(name, "dCOR - |PCC|"))
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
  R$COV <- label.dep.m(cov(t(M)), labels)

  # dcor hclust order
  D.DCOR.r <- as.dist(1-cor(t(R$DCOR), method="pearson")) + dist(R$DCOR)*2
  Rowv <- rowMeans(R$DCOR, na.rm = TRUE)
  Rhclust <- as.dendrogram(hclust(D.DCOR.r, method="average"))
  Rhclust <- reorder(Rhclust, Rowv)
  R$Rhclust
  R$rowInd <- order.dendrogram(Rhclust)
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

plot.sploms <- function(CLS, D, name, MIN=0.18) {
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
  dev.off()

  pdf(paste0("splom.summary.",name,".pdf"), width=12, height=12)
  summary.plots(CLS$CLS,D$DCOR, sym=T)
  dev.off()
  CD.rowInd
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
  plot.dep.matrix(R$GEO.D, name=paste0("GEO.GOLD.4mut",name))
  
  plot.pairs(M=exprs(SCAN[gold.i,]), R=R$SCAN.D, labels=gold.labels, name=paste0("SCAN.GOLD.",name))
  plot.pairs(M=exprs(GEO[gold.i,]), R=R$GEO.D, labels=gold.labels, name=paste0("GEO.GOLD.",name))
  
  R$SCAN.CLS <- calc.plot.classes(M=exprs(SCAN[gold.i,]), b=GSE2180.SCAN.b, labels=gold.labels, name=paste0("SCAN.GOLD.",name))
  R$GEO.CLS <- calc.plot.classes(M=exprs(GEO[gold.i,]), b=GSE2180.GEO.b, labels=gold.labels, name=paste0("GEO.GOLD.",name))
  
  R$SCAN.splom.rowInd <- plot.sploms(CLS=R$SCAN.CLS, D=R$SCAN.D, name=paste0("SCAN.GOLD.",name))
  R$GEO.splom.rowInd <- plot.sploms(CLS=R$GEO.CLS, D=R$GEO.D, name=paste0("GEO.GOLD.",name))
  R
}
## ==================================================



# 4mut
all.4mut <- analyze.and.plot(GEO=GSE2180.GEO, SCAN=GSE2180.SCAN, name="4mut")

# WT+ms
qq <- GSE2180.GEO$genotype %in% c('N2','ms')
wt.ms.2mut <- analyze.and.plot(GEO=GSE2180.GEO[,qq], SCAN=GSE2180.SCAN[,qq], name="WT+ms")

