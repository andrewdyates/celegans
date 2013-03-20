library(sva)
library(lumi) # for nice figures, plus Biobase
library(Biobase)
library(energy)
library("RColorBrewer")
source("~/Dropbox/biostat/git_repos/dependency_glyph_splom/lib.R")
source("~/Dropbox/biostat/git_repos/boolean_implication_fit/bool.R")
source("~/Dropbox/biostat/git_repos/boolean_implication_fit/step.up.R")

all.pairs.dcor <- function(M) {
  n <- nrow(M)
  D <- outer(1:n,1:n, FUN = Vectorize(function(i,j) dcor(M[i,],M[j,])))
  rownames(D) <- rownames(M)
  colnames(D) <- rownames(M)
  D
}

all.steps <- function(M, do.plot=TRUE) {
  STEPS <- apply(M, 1, fit.upstep)
  if(do.plot) {
    for (i in 1:dim(M)[1]) {
      title <- rownames(M)[i], names(unlist(row.nums))[i])
      plot.stepfit(STEPS[[i]], v=M[i,], add.mean.median=T, main=paste(title, "Stepfit"))
      plot.sse(STEPS[[i]], add.mean.median=T, main=paste(title, "SSE"))
    }
  }
  STEPS
}

all.pairs.cls <- function(M, steps, do.plot=TRUE) {
  n <- dim(M)[1]
  CLS <- mat.or.vec(n,n)
  for(i in 1:n) { # row
    for(j in 1:n) { # col
      y <- M[i,]
      y.th <- steps[[i]]$th
      x <- M[j,]
      x.th <- steps[[j]]$th
      x.title=paste(rownames(M)[j], names(unlist(row.nums))[j])
      y.title=paste(rownames(M)[i], names(unlist(row.nums))[i])
      RR <- cls.pair(x,y,x.th,y.th, b.x=b, b.y=b, do.plot=(i<j)&&do.plot, xlab=x.title, ylab=y.title)
      CLS[i,j] <- cls.to.enum(RR$CLS)
    }
  }
  rownames(CLS) <- rownames(M)
  colnames(CLS) <- rownames(M)
  CLS
}


gold.genes <- c("pal-1", "tbx-8", "tbx-9", "elt-1", "cwn-1", "hnd-1", "scrt-1", "elt-3", "nhr-25", "nob-1", "vab-7", "unc-120", "hlh-1", "mab-21", "lin-26")

load("../GSE2180.ALL.NOFILT.RData")
load("../GSE9665.ALL.NOFILT.RData")
load("../geo.GSE2180.GSE9665.RData")
GPL200 <- read.table("../GPL200-2880.txt", sep="\t", quote="", comment="", header=TRUE, row.names=1)

## ========================================
## Row align scan normalized to GEO GPL row order

all(rownames(GSE2180.ALL.SCAN) == rownames(GSE2180.ALL.UPC))
qq <- match(rownames(exprs(GSE2180)), rownames(exprs(GSE2180.ALL.SCAN)))
all(qq == 1:dim(GSE2180)[1])
# [1] FALSE
any(is.na(qq))
#[1] FALSE
GSE2180.ALL.SCAN <- GSE2180.ALL.SCAN[qq,]
GSE2180.ALL.UPC <- GSE2180.ALL.UPC[qq,]

all(rownames(exprs(GSE9665.ALL.SCAN)) == rownames(exprs(GSE9665.ALL.UPC)))
qq <- match(rownames(exprs(GSE9665)), rownames(exprs(GSE9665.ALL.SCAN)))
all(qq == 1:dim(GSE9665)[1])
# [1] FALSE
any(is.na(qq))
#[1] FALSE
GSE9665.ALL.SCAN <- GSE9665.ALL.SCAN[qq,]
GSE9665.ALL.UPC <- GSE9665.ALL.UPC[qq,]

all(rownames(exprs(GSE9665)) == rownames(exprs(GSE2180)))
# [1] TRUE
all(rownames(exprs(GSE9665)) == rownames(exprs(GSE2180.ALL.SCAN)))
# [1] TRUE
## OK: row IDs align

## Column IDs (samples) must also align
colnames(exprs(GSE2180.ALL.SCAN)) <- sub(".CEL.gz", "", colnames(exprs(GSE2180.ALL.SCAN)))
colnames(exprs(GSE2180.ALL.UPC)) <- sub(".CEL.gz", "", colnames(exprs(GSE2180.ALL.UPC)))
colnames(exprs(GSE9665.ALL.SCAN)) <- sub(".cel.gz", "", colnames(exprs(GSE9665.ALL.SCAN)))
colnames(exprs(GSE9665.ALL.UPC)) <- sub(".cel.gz", "", colnames(exprs(GSE9665.ALL.UPC)))

qq <- match(colnames(exprs(GSE2180)), colnames(exprs(GSE2180.ALL.SCAN)))
all(qq == 1:dim(GSE2180)[2])
qq <- match(colnames(exprs(GSE2180)), colnames(exprs(GSE2180.ALL.UPC)))
all(qq == 1:dim(GSE2180)[2])
qq <- match(colnames(exprs(GSE9665)), colnames(exprs(GSE9665.ALL.SCAN)))
all(qq == 1:dim(GSE9665)[2])
qq <- match(colnames(exprs(GSE9665)), colnames(exprs(GSE9665.ALL.UPC)))
all(qq == 1:dim(GSE9665)[2])
## OK: All aligned.

## Filter possibly low qualtiy probes from GSE2180
low.qual.arrays <- c("GSM39507")
col.filt <- !colnames(exprs(GSE2180)) %in% low.qual.arrays
GSE2180 <- GSE2180[,col.filt]
GSE2180.ALL.SCAN <- GSE2180.ALL.SCAN[,col.filt]
GSE2180.ALL.UPC <- GSE2180.ALL.UPC[,col.filt]
## [1] 

## ----------------------------------------
## GOLD STANDARD NETWORKS
## ----------------------------------------
# GSE2180 Time Series
row.nums <- sapply(gold.genes, function(s) grep(paste0("(",s," |",s,"$)"), featureData(GSE2180)$Gene.Symbol))
rows <- unlist(row.nums)
## visually inspect, get correct probe IDs
sapply(row.nums, function(q) featureData(GSE2180)$Gene.Symbol[q]) ## OK
GSE2180.GEO.GOLD <- GSE2180[unlist(row.nums),]
GSE2180.SCAN.GOLD <- GSE2180.ALL.SCAN[unlist(row.nums),]

## ----------------------------------------
## DEPENDENCIES
## ----------------------------------------
GSE2180.GEO.GOLD.DCOR <- all.pairs.dcor(exprs(GSE2180.GEO.GOLD))
GSE2180.GEO.GOLD.PCC <- cor(t(exprs(GSE2180.GEO.GOLD)))

GSE2180.SCAN.GOLD.DCOR <- all.pairs.dcor(exprs(GSE2180.SCAN.GOLD))
GSE2180.SCAN.GOLD.PCC <- cor(t(exprs(GSE2180.SCAN.GOLD)))


## GET BOOL CLASSES
## ------------------------------
pdf("step.plots.pdf")
STEPS <- all.steps(exprs(GSE2180.SCAN.GOLD))
dev.off()

# remove background probes
upc.maxs <- apply(exprs(GSE2180.ALL.UPC),1,max)
stds <- apply(exprs(GSE2180.ALL.SCAN[upc.maxs>0.3,]),1,sd)
b <- quantile(stds,0.03)*2
## b == 0.1726519 

pdf(paste0("bool.scatterplots.all.pdf"))
GSE2180.SCAN.GOLD.CLS <- all.pairs.cls(GSE2180.ALL.SCAN, steps=STEPS)
dev.off()

## ## TEST SINGLE
## i=1; j=4
## y <- M[i,]
## y.th <- STEPS[[i]]$th
## x <- M[j,]
## x.th <- STEPS[[j]]$th
## x.title=paste(rownames(M)[j], names(unlist(row.nums))[j])
## y.title=paste(rownames(M)[i], names(unlist(row.nums))[i])
## pdf("test.pdf"); cls.pair(x,y,x.th,y.th, b.x=b, b.y=b, do.plot=i<j, xlab=x.title, ylab=y.title); dev.off()

# ========================================

# histograms
pdf("dCOR.PCC.GEO.SCAN.hist.pdf")
hist(GSE2180.SCAN.GOLD.DCOR-abs(GSE2180.SCAN.GOLD.PCC), main="SCAN dCOR - SCAN |PCC|")
hist(GSE2180.SCAN.GOLD.DCOR-GSE2180.GEO.GOLD.DCOR, main="SCAN dCOR - GEO dCOR")
hist(GSE2180.GEO.GOLD.DCOR-abs(GSE2180.GEO.GOLD.PCC), main="GEO dCOR - GEO |PCC|")
hist(GSE2180.SCAN.GOLD.DCOR-abs(GSE2180.GEO.GOLD.PCC), main="SCAN dCOR - GEO |PCC|")
hist(abs(GSE2180.SCAN.GOLD.PCC)-abs(GSE2180.GEO.GOLD.PCC), main="SCAN |PCC| - GEO |PCC|")
dev.off()

MIN <- 0
MAX <- 1
cols <- brewer.pal(8,"RdYlBu")
heatmap_breaks <- seq(MIN,MAX,0.01)
heatmap_cols <- rev(colorRampPalette(cols)(length(heatmap_breaks)-1))

heatmap_breaks.2 <- seq(-0.4,0.4,0.01)
heatmap_cols.2 <- rev(colorRampPalette(cols)(length(heatmap_breaks.2)-1))

pdf("DCOR.ALL.pdf", width=20, height=20)
R <- heatmap.3(GSE2180.SCAN.GOLD.DCOR, MIN=0, MAX=1, symm=T, main="SCAN DCOR")
heatmap.2(abs(GSE2180.SCAN.GOLD.PCC), main="SCAN |PCC|", symm=T, col=heatmap_cols, breaks=heatmap_breaks, Colv=R$colDendrogram, Rowv=R$rowDendrogram, trace='none')
heatmap.2(abs(GSE2180.GEO.GOLD.DCOR), main="GEO DCOR", symm=T, col=heatmap_cols, breaks=heatmap_breaks, Colv=R$colDendrogram, Rowv=R$rowDendrogram, trace='none')
heatmap.2(abs(GSE2180.GEO.GOLD.PCC), main="GEO |PCC|", symm=T, col=heatmap_cols, breaks=heatmap_breaks, Colv=R$colDendrogram, Rowv=R$rowDendrogram, trace='none')

heatmap.2(GSE2180.SCAN.GOLD.DCOR - abs(GSE2180.SCAN.GOLD.PCC), main="SCAN DCOR - SCAN |PCC|", symm=T, col=heatmap_cols.2, breaks=heatmap_breaks.2, Colv=R$colDendrogram, Rowv=R$rowDendrogram, trace='none')
heatmap.2(GSE2180.SCAN.GOLD.DCOR - abs(GSE2180.GEO.GOLD.PCC), main="SCAN DCOR - GEO |PCC|", symm=T, col=heatmap_cols.2, breaks=heatmap_breaks.2, Colv=R$colDendrogram, Rowv=R$rowDendrogram, trace='none')
heatmap.2(abs(GSE2180.SCAN.GOLD.PCC) - abs(GSE2180.GEO.GOLD.PCC), main="SCAN |PCC| - GEO |PCC|", symm=T, col=heatmap_cols.2, breaks=heatmap_breaks.2, Colv=R$colDendrogram, Rowv=R$rowDendrogram, trace='none')
pairs(t(exprs(GSE2180.SCAN.GOLD[R$rowInd,])), main="SCAN GOLD")
pairs(t(exprs(GSE2180.SCAN.GOLD[R$rowInd,])), labels=names(unlist(row.nums))[R$rowInd])
dev.off()



## ----------------------------------------
## PLOTS
## ----------------------------------------
# Plot differences in normalization
pdf("gse2180.gold.geo.vs.scan.plots.pdf")
for (i in 1:length(rows)) {
  x <- exprs(GSE2180.GEO.GOLD)[i,]
  y <- exprs(GSE2180.SCAN.GOLD)[i,]
  rho <- cor(x, y)
  title <- paste(names(unlist(row.nums))[i], rownames(exprs(GSE2180.GEO.GOLD))[i], format(rho,2))
  plot(x,y,main=title, xlab="GEO RMA norm", ylab="SCAN norm")
}
dev.off()

# Plot each gene in time and for mutants per norm procedure
time <- GSE2180$time
pdf("gse2180.time.mutant.gold.geo.vs.scan.plots.pdf")
i <- 1
for (i in 1:length(rows)) {
  x <- exprs(GSE2180.GEO.GOLD)[i,]
  y <- exprs(GSE2180.SCAN.GOLD)[i,]
  rho <- cor(x, y)
  title <- paste(names(unlist(row.nums))[i], rownames(exprs(GSE2180.GEO.GOLD))[i], "rho=", format(rho,digits=4))
  par(mar=c(5,4,4,5)+.1)
  plot(time, x, col="#cc0000", pch=as.numeric(GSE2180$genotype), xlab="time", ylab="SCAN Intensity", main=title)
  par(new=TRUE)
  plot(time, y, col="#000099", xaxt="n",yaxt="n",xlab="",ylab="",pch=as.numeric(GSE2180$genotype))
  axis(4)
  mtext("RMA Intensity",side=4,line=3)
  legend("topleft",col=c("#cc0000","#000099"),lty=1,legend=c("RMA","SCAN"))
  legend("topright",pch=1:4,legend=levels(GSE2180$genotype))
}
dev.off()


#> levels(GSE2180$genotype)
#[1] "ms" "N2" "pi" "pp"
# Plot each mutant in time per gene
time <- GSE2180$time
col.scale <- brewer.pal(4,"Set1")
cols <- col.scale[as.numeric(GSE2180$genotype)]
pdf("gse2180.time.mutant.gold.geo.vs.scan.plots.pdf")
for (i in 1:length(rows)) {
  y <- exprs(GSE2180.SCAN.GOLD)[i,]
  title <- paste(names(unlist(row.nums))[i], rownames(exprs(GSE2180.GEO.GOLD))[i])
  plot(time, y, col=cols, pch=20, main=title, ylab="SCAN Intensity")
  legend("topleft", pch=20, col=col.scale, legend=levels(GSE2180$genotype))
}
dev.off()


pdf("gse2180.geo.gold.pairs.pdf", width=30, height=30)
pairs(t(exprs(GSE2180.GEO.GOLD)), labels=names(unlist(row.nums)), main="GEO GSE2180 C.Elegans GSN Gene Names (with multi-probe enumerations)")
pairs(t(exprs(GSE2180.GEO.GOLD)), main="GEO GSE2180 C.Elegans GSN Probes")
pairs(t(exprs(GSE2180.SCAN.GOLD)), labels=names(unlist(row.nums)), main="SCAN GSE2180 C.Elegans GSN Gene Names (with multi-probe enumerations)")
pairs(t(exprs(GSE2180.SCAN.GOLD)), main="SCAN GSE2180 C.Elegans GSN Probes")
dev.off()

# ----------
# GSE9665 Knockouts
GSE9665.GEO.GOLD <- GSE9665[unlist(row.nums),]
GSE9665.SCAN.GOLD <- GSE9665.ALL.SCAN[unlist(row.nums),]

# Plot differences in normalization
pdf("gse9665.gold.geo.vs.scan.plots.pdf")
for (i in 1:length(rows)) {
  x <- exprs(GSE9665.GEO.GOLD)[i,]
  y <- exprs(GSE9665.SCAN.GOLD)[i,]
  rho <- cor(x, y)
  title <- paste(names(unlist(row.nums))[i], rownames(exprs(GSE9665.GEO.GOLD))[i], format(rho,2))
  plot(x,y,main=title, xlab="GEO RMA norm", ylab="SCAN norm")
}
dev.off()

pdf("gse9665.geo.scan.gold.pairs.pdf", width=30, height=30)
  pairs(t(exprs(GSE9665.GEO.GOLD)), labels=names(unlist(row.nums)), main="RMA GSE9665 C.Elegans GSN Gene Names (with multi-probe enumerations)")
  pairs(t(exprs(GSE9665.GEO.GOLD)), main="RMA GSE9665 C.Elegans GSN Probes")
  pairs(t(exprs(GSE9665.SCAN.GOLD)), labels=names(unlist(row.nums)), main="SCAN GSE9665 C.Elegans GSN Gene Names (with multi-probe enumerations)")
  pairs(t(exprs(GSE9665.SCAN.GOLD)), main="SCAN GSE9665 C.Elegans GSN Probes")
dev.off()
## ========================================



## GLYPH SPLOM
## ----------------------------------------

## manually order
#GSE2180.SCAN.GOLD.CLS, GSE2180.SCAN.GOLD.DCOR



pdf("testsplom.pdf")
splom(GSE2180.SCAN.GOLD.CLS, GSE2180.SCAN.GOLD.DCOR, asGlyphs=T, MAX=1, MIN=0.25, reorder=F)
dev.off()

pdf("test.summary.pdf")
summary.plots(GSE2180.SCAN.GOLD.CLS,GSE2180.SCAN.GOLD.DCOR, sym=T)
dev.off()


M <- exprs(GSE2180.SCAN.GOLD)
CLS <- GSE2180.SCAN.GOLD.CLS
DCOR <- GSE2180.SCAN.GOLD.DCOR
D.DCOR.r <- as.dist(1-cor(t(DCOR), method="pearson")) + dist(DCOR)*2
D.cls.r <- dist(CLS)
Rowv <- rowMeans(DCOR, na.rm = TRUE)

DCOR.weight=2
Rhclust <- as.dendrogram(hclust(D.DCOR.r*DCOR.weight+D.cls.r, method="average"))
Rhclust <- reorder(Rhclust, Rowv)
rowInd <- order.dendrogram(Rhclust)

pdf("raghu.mar12.gsplom.pdf", width=20,height=20)
splom(CLS[rowInd,rowInd], DCOR[rowInd,rowInd], asGlyphs=T, MAX=1, MIN=0.17, reorder=F)
pairs(t(M[rowInd,]))
pairs(t(M[rowInd,]),labels=names(unlist(row.nums))[rowInd])
dev.off()


## Filter genotype for WT, mm
## ----------------------------------------
all(match(sub(".CEL.gz","",sampleNames(GSE2180.SCAN.GOLD)),sampleNames(GSE2180.GEO.GOLD)) == 1:122)
all(match(sub(".CEL.gz","",colnames(exprs(GSE2180.SCAN.GOLD))),sampleNames(GSE2180.GEO.GOLD)) == 1:122)


C <- GSE2180.SCAN.GOLD[,GSE2180.GEO.GOLD$genotype %in% c('ms','N2')]
C.DCOR <- all.pairs.dcor(exprs(C))
C.PCC <- cor(t(exprs(C)))

pdf("step.plots.C.pdf")
C.STEPS <- all.steps(exprs(C))
dev.off()

pdf("bool.scatterplots.C.pdf")
C.CLS <- all.pairs.cls(exprs(C),steps=C.STEPS)
dev.off()

pdf("C.splom.summary.pdf")
summary.plots(C.CLS,C.DCOR, sym=T)
dev.off()

D.DCOR.r <- as.dist(1-cor(t(C.DCOR), method="pearson")) + dist(C.DCOR)*2
D.cls.r <- dist(C.CLS)
Rowv <- rowMeans(C.DCOR, na.rm = TRUE)
DCOR.weight=2
Rhclust <- as.dendrogram(hclust(D.DCOR.r*DCOR.weight+D.cls.r, method="average"))
Rhclust <- reorder(Rhclust, Rowv)
rowInd <- order.dendrogram(Rhclust)

pdf("C.mar14.gsplom.pdf", width=20,height=20)
splom(C.CLS[rowInd,rowInd], C.DCOR[rowInd,rowInd], asGlyphs=T, MAX=1, MIN=0.17, reorder=F)
splom(C.CLS[rowInd,rowInd], C.DCOR[rowInd,rowInd], asGlyphs=T, MAX=1, MIN=0.3, reorder=F)
pairs(t(exprs(C)[rowInd,]))
pairs(t(exprs(C)[rowInd,]),labels=names(unlist(row.nums))[rowInd])
dev.off()


pdf("dCOR.PCC.C.SCAN.hist.pdf")
hist(C.DCOR, main="SCAN dCOR C")
hist(C.PCC, main="SCAN PCC C")
hist(C.DCOR - GSE2180.SCAN.GOLD.DCOR, main="SCAN dCOR C - SCAN dCOR All")
hist(abs(C.PCC) - abs(GSE2180.SCAN.GOLD.PCC), main="SCAN |PCC| C - SCAN |PCC| All")
hist(C.DCOR - abs(C.PCC), main="SCAN dCOR C - SCAN |PCC| C")
dev.off()

pdf("DCOR.C.pdf", width=20, height=20)
R <- heatmap.3(C.DCOR, MIN=0, MAX=1, symm=T, main="SCAN DCOR C")
heatmap.2(C.PCC, main="SCAN PCC C", symm=T, col=heatmap_cols, breaks=heatmap_breaks, Colv=R$colDendrogram, Rowv=R$rowDendrogram, trace='none')
heatmap.2(abs(C.PCC), main="SCAN |PCC| C", symm=T, col=heatmap_cols, breaks=heatmap_breaks, Colv=R$colDendrogram, Rowv=R$rowDendrogram, trace='none')

heatmap.2(C.DCOR - GSE2180.SCAN.GOLD.DCOR, main="SCAN DCOR C - SCAN DCOR ALL", symm=T, col=heatmap_cols.2, breaks=heatmap_breaks.2, Colv=R$colDendrogram, Rowv=R$rowDendrogram, trace='none')
heatmap.2(C.DCOR - abs(C.PCC), main="SCAN DCOR C - SCAN |PCC| ALL", symm=T, col=heatmap_cols.2, breaks=heatmap_breaks.2, Colv=R$colDendrogram, Rowv=R$rowDendrogram, trace='none')
heatmap.2(abs(C.PCC) - abs(GSE2180.SCAN.GOLD.PCC), main="SCAN |PCC| C - SCAN |PCC| ALL", symm=T, col=heatmap_cols.2, breaks=heatmap_breaks.2, Colv=R$colDendrogram, Rowv=R$rowDendrogram, trace='none')
pairs(t(exprs(C[R$rowInd,])), main="SCAN GOLD C probes")
pairs(t(exprs(C[R$rowInd,])), main="SCAN GOLD C genes", labels=names(unlist(row.nums))[R$rowInd])
dev.off()




# ========================================
C <- GSE2180.SCAN.GOLD[,GSE2180.GEO.GOLD$genotype %in% c('ms')]
C.DCOR <- all.pairs.dcor(exprs(C))
C.PCC <- cor(t(exprs(C)))

pdf("step.plots.Cms.pdf")
C.STEPS <- all.steps(exprs(C))
dev.off()

pdf("bool.scatterplots.Cms.pdf")
C.CLS <- all.pairs.cls(exprs(C),steps=C.STEPS)
dev.off()

pdf("Cms.splom.summary.pdf")
summary.plots(C.CLS,C.DCOR, sym=T)
dev.off()

D.DCOR.r <- as.dist(1-cor(t(C.DCOR), method="pearson")) + dist(C.DCOR)*2
D.cls.r <- dist(C.CLS)
Rowv <- rowMeans(C.DCOR, na.rm = TRUE)
DCOR.weight=2
Rhclust <- as.dendrogram(hclust(D.DCOR.r*DCOR.weight+D.cls.r, method="average"))
Rhclust <- reorder(Rhclust, Rowv)
rowInd <- order.dendrogram(Rhclust)

pdf("Cms.mar14.gsplom.pdf", width=20,height=20)
splom(C.CLS[rowInd,rowInd], C.DCOR[rowInd,rowInd], asGlyphs=T, MAX=1, MIN=0.17, reorder=F)
splom(C.CLS[rowInd,rowInd], C.DCOR[rowInd,rowInd], asGlyphs=T, MAX=1, MIN=0.3, reorder=F)
pairs(t(exprs(C)[rowInd,]))
pairs(t(exprs(C)[rowInd,]),labels=names(unlist(row.nums))[rowInd])
dev.off()


pdf("dCOR.PCC.Cms.SCAN.hist.pdf")
hist(C.DCOR, main="SCAN dCOR Cms")
hist(C.PCC, main="SCAN PCC Cms")
hist(C.DCOR - GSE2180.SCAN.GOLD.DCOR, main="SCAN dCOR Cms - SCAN dCOR All")
hist(abs(C.PCC) - abs(GSE2180.SCAN.GOLD.PCC), main="SCAN |PCC| Cms - SCAN |PCC| All")
hist(C.DCOR - abs(C.PCC), main="SCAN dCOR Cms - SCAN |PCC| Cms")
dev.off()

pdf("DCOR.Cms.pdf", width=20, height=20)
R <- heatmap.3(C.DCOR, MIN=0, MAX=1, symm=T, main="SCAN DCOR Cms")
heatmap.2(C.PCC, main="SCAN PCC Cms", symm=T, col=heatmap_cols, breaks=heatmap_breaks, Colv=R$colDendrogram, Rowv=R$rowDendrogram, trace='none')
heatmap.2(abs(C.PCC), main="SCAN |PCC| Cms", symm=T, col=heatmap_cols, breaks=heatmap_breaks, Colv=R$colDendrogram, Rowv=R$rowDendrogram, trace='none')

heatmap.2(C.DCOR - GSE2180.SCAN.GOLD.DCOR, main="SCAN DCOR Cms - SCAN DCOR ALL", symm=T, col=heatmap_cols.2, breaks=heatmap_breaks.2, Colv=R$colDendrogram, Rowv=R$rowDendrogram, trace='none')
heatmap.2(C.DCOR - abs(C.PCC), main="SCAN DCOR Cms - SCAN |PCC| ALL", symm=T, col=heatmap_cols.2, breaks=heatmap_breaks.2, Colv=R$colDendrogram, Rowv=R$rowDendrogram, trace='none')
heatmap.2(abs(C.PCC) - abs(GSE2180.SCAN.GOLD.PCC), main="SCAN |PCC| Cms - SCAN |PCC| ALL", symm=T, col=heatmap_cols.2, breaks=heatmap_breaks.2, Colv=R$colDendrogram, Rowv=R$rowDendrogram, trace='none')
pairs(t(exprs(C[R$rowInd,])), main="SCAN GOLD Cms probes")
pairs(t(exprs(C[R$rowInd,])), main="SCAN GOLD Cms genes", labels=names(unlist(row.nums))[R$rowInd])
dev.off()
