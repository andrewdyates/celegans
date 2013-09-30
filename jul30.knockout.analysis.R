Tab.network.genes <- read.table("../c.elegans.knockdown.edges.network.genes.csv", sep=",", header=T, stringsAsFactors=F)
Tab.ko.edges <- read.table("../c.elegans.knockdown.edges.edge.table.csv", sep=",", header=T, stringsAsFactors=F)

Activators <- Tab.ko.edges[Tab.ko.edges$edgetype==1,]
Repressors <- Tab.ko.edges[Tab.ko.edges$edgetype==-1,]

# NOTES: tbx-8/9 knocked out together, vab-7 not detected, lin-26 discarded due to unreliable probe
sort(summary(as.factor(Activators$knockout)), decreasing=T)
#  pal-1   tbx-8   tbx-9   elt-1   hnd-1  scrt-1   hlh-1  nhr-25   nob-1   elt-3 
#      8       7       7       4       4       4       3       3       2       1 
# lin-26 unc-120   vab-7 
#      1       1       1
# AUTO-EDGES
self.act <- Activators$knockout==Activators$genename
Activators$knockout[self.act]
#[1] "tbx-8" "elt-1" "hnd-1" "pal-1"

# ALL EDGES IN NETWORK
# A:pal-1, B:tbx-8, B:hnd-1, C:tbx-9, E:elt-1, F:unc-120, G:scrt-1, H:nob-1, I:elt-3, J:nhr-25, J:hlh-1

sort(summary(as.factor(Repressors$knockout)), decreasing=T)
# nhr-25   pal-1   nob-1   elt-3   vab-7   tbx-8   tbx-9   hlh-1 unc-120  lin-26 
#     54      54      43      33      27      25      25      23      22      21 
# scrt-1   elt-1   hnd-1 
#     15      13      10

# WHICH GENES CONTINUE TO REAPPEAR? (note: 13 knockouts total)
sort(summary(as.factor(Activators$genename)), decreasing=T)
       ##  tbx-8       nhr-188 CELE_F23B12.7         cwn-1         hnd-1 
       ##     13             4             3             3             2 
       ##  mes-2       nhr-104       nhr-163       nhr-171       nhr-205 
       ##      2             2             2             2             2 
       ## nhr-86  CELE_T13F2.2         elt-1       F58G1.2         hlh-1 
       ##      2             1             1             1             1 
       ##  lim-4       nhr-222        nhr-61         nob-1         pal-1 
       ##      1             1             1             1             1 
sort(summary(as.factor(Activators$genename[!self.act])), decreasing=T)
 ##        tbx-8       nhr-188 CELE_F23B12.7         cwn-1         mes-2 
 ##           12             4             3             3             2 
 ##      nhr-104       nhr-163       nhr-171       nhr-205        nhr-86 
 ##            2             2             2             2             2 
 ## CELE_T13F2.2       F58G1.2         hlh-1         hnd-1         lim-4 
 ##            1             1             1             1             1 
 ##      nhr-222        nhr-61         nob-1 
 ##            1             1             1 
sort(summary(as.factor(Activators$genename[!self.act&Activators$knock!="tbx-9"])), decreasing=T)
      ##   tbx-8       nhr-188 CELE_F23B12.7         cwn-1       nhr-171 
      ##      11             4             2             2             2 
      ## nhr-205  CELE_T13F2.2       F58G1.2         hlh-1         hnd-1 
      ##       2             1             1             1             1 
      ##   lim-4         mes-2       nhr-104       nhr-163       nhr-222 
      ##       1             1             1             1             1 
      ##  nhr-61        nhr-86         nob-1 
      ##       1             1             1 
sort(summary(as.factor(Activators$genename[!self.act&Activators$knock!="tbx-9"&!is.na(Activators$bool)])), decreasing=T)
      ##   tbx-8 CELE_F23B12.7         cwn-1       nhr-171  CELE_T13F2.2 
      ##       9             2             2             2             1 
      ## F58G1.2         hlh-1         hnd-1         mes-2       nhr-104 
      ##       1             1             1             1             1 
      ## nhr-163       nhr-222        nhr-61        nhr-86         nob-1 
      ##       1             1             1             1             1


## BARPLOT FOR BOOL
qq <- !self.act&Activators$knock!="tbx-9"&!is.na(Activators$bool)
barplot(summary(as.factor(Activators[qq,]$bool)))
barplot(summary(as.factor(Activators[qq,]$weak)))
dim(BOOL.TF.F) # 655
n <- sum(qq) # 26
bool.pct.act <- summary(as.factor(Activators[qq,]$bool))/n*100
bool.pct.act[["6"]] <- 0
bool.pct.act[["7"]] <- 0
qq <- Repressors$knock!="tbx-9"&!is.na(Repressors$bool)
n <- sum(qq) # 245
bool.pct.rep <- summary(as.factor(Repressors[qq,]$bool))/n*100
bool.pct.rep[["6"]] <- 0
bool.pct.rep[["7"]] <- 0
n <- sum(upper.tri(BOOL.TF.F))*2 # 428370
bool.pct.all <- summary(as.factor(BOOL.TF.F[!diag(ncol(BOOL.TF.F))]))/n*100
# make this pretty
pdf("~/Desktop/boolclass.barchart.pdf", width=8, height=8)
barplot(rbind(bool.pct.act,bool.pct.rep,bool.pct.all), beside=T, xlab="Abundance Level Boolean Class: Knockout (KO) on X-axis versus Target (Y) on Y-axis", ylab="Percent Total", names.arg=c("1: KO<=Y","2: KO<=>Y","3: KO=>Y","4: UNL","5: MX","6: NC","7: OR"), legend.text=c("Activating (26)", "Repressing (245)", "All (428k)"))
dev.off()

#BARPLOT FOR WEAK
qq <- !self.act&Activators$knock!="tbx-9"&!is.na(Activators$weak)
n <- sum(qq) # 26
weak.pct.act <- summary(as.factor(Activators[qq,]$weak))/n*100
weak.pct.act[["1"]] <- 0
weak.pct.act[["4"]] <- 0
qq <- Repressors$knock!="tbx-9"&!is.na(Repressors$weak)
n <- sum(qq) # 245
weak.pct.rep <- summary(as.factor(Repressors[qq,]$weak))/n*100
n <- sum(upper.tri(BOOL.TF.F))*2 # 428370
weak.pct.all <- summary(as.factor(WEAK.TF.F[!diag(ncol(WEAK.TF.F))]))/n*100
# make this pretty
pdf("~/Desktop/weakclass.barchart.pdf", width=8, height=8)
barplot(rbind(weak.pct.act,weak.pct.rep,weak.pct.all), beside=T, xlab="Detection Level Class: Knockout (KO) versus Target (Y)", ylab="Percent Total", names.arg=c("1: AND","2: KO<=Y","3: KO=>Y","4: XOR","Other"), legend.text=c("Activating (26)", "Repressing (245)", "All (428k)"))
dev.off()

#BOXPLOT FOR dCor
qq <- !self.act&Activators$knock!="tbx-9"&!is.na(Activators$dcor)
n <- sum(qq) # 26
dcor.act <- Activators[qq,]$dcor
qq <- Repressors$knock!="tbx-9"&!is.na(Repressors$dcor)
n <- sum(qq) # 245
dcor.rep <- Repressors[qq,]$dcor
n <- sum(upper.tri(BOOL.TF.F))*2 # 428370
dcor.all <- DCOR.TF.F[!diag(ncol(DCOR.TF.F))]
# make this pretty
pdf("~/Desktop/dcor.boxplot.pdf", width=8, height=8)
boxplot(dcor.act,dcor.rep,dcor.all,names=c("Activating (26)", "Repressing (245)", "All (428k)"), ylab="dCor")
dev.off()


Activators[!self.act&Activators$knock!="tbx-9"&!is.na(Activators$bool),]
## 16    scrt-1 176450        1         tbx-8    4    2 0.6740  
## 17    scrt-1 179947        1 CELE_F23B12.7    4    5 0.2572  unknown lower cluster, undefined relation
## 45     tbx-8 175096        1         mes-2    5    5 0.6705  pie-1, pos-1 cluster
## 47     tbx-8 177248        1       nhr-104    4    5 0.2003  top "all blue" cluster that is inverse to dcp−66
## 48     tbx-8 178685        1        nhr-86    4    5 0.3610  near CELE_F23B12.7, but loosely related
## 49     tbx-8 179947        1 CELE_F23B12.7    4    5 0.1853  
## 50     tbx-8 183182        1       nhr-163    3    3 0.6178  mab-21 cluster
## 51     tbx-8 173399        1         cwn-1    4    2 0.7631
## 97     elt-1 175142        1        nhr-61    4    5 0.3823  unknown lower cluster, unrelated to phase 2, but necessary for phase 3 
## 98     elt-1 176450        1         tbx-8    2    3 0.9235
## 100    elt-1 183807        1       nhr-171    3    3 0.5664  near phase-3 cluster
## 124    hlh-1 176450        1         tbx-8    1    2 0.6843
## 126    hlh-1 188475        1  CELE_T13F2.2    5    2 0.4813  pos-1, pie-1 cluster
## 181   nhr-25 176450        1         tbx-8    1    2 0.7270
## 217    elt-3 176450        1         tbx-8    1    2 0.6744
## 261    nob-1 174931        1       F58G1.2    5    5 0.4630  unknown cluster that depends on pos-1, pie-1 cluster
## 262    nob-1 176450        1         tbx-8    1    2 0.7467
## 285  unc-120 176450        1         tbx-8    1    3 0.7984
## 296    hnd-1 176450        1         tbx-8    2    3 0.8993
## 354    pal-1 173788        1         hlh-1    4    3 0.2015
## 356    pal-1 176450        1         tbx-8    4    3 0.2085
## 357    pal-1 176641        1         nob-1    4    3 0.2016
## 358    pal-1 183457        1         hnd-1    4    3 0.2229
## 359    pal-1 183807        1       nhr-171    4    3 0.3286
## 360    pal-1 188828        1       nhr-222    4    3 0.2944  top "all blue" cluster that is inverse to dcp−66
## 361    pal-1 173399        1         cwn-1    4    3 0.4244

## ----------------------------------------
# are there any nodes in the network in the top edge list?
unique(Tab.ko.edges$genename[Tab.ko.edges$genename %in% Tab.network.genes$gene])
#[1] "tbx-8"        "cwn-1"        "elt-1"        "CELE_M03D4.4" "ldb-1"       
#[6] "hnd-1"        "hlh-1"        "pal-1"        "nob-1"       
# CELE_M03D4.4 (mab-21 K cluster)
# ldb-1 (elt-1 E cluster)

# ------------------------------
length(unique(Repressors$genename))
# [1] 126 of ~700, 2 in 19 of unknown genes
#ldb-1: 2 # nob-1, nhr-25
#CELE_M03D4.4: 2 # nob-1, nhr-25
## lin-32: top repressor (11), top part of gold clust, nec on all but not related to mab-21
## chd-1: pos-1 clust
## K09A11.1: unknown bottom clust

# -----------------------------
summary(as.factor(Activators$bool[!self.act]))
# histogram of bool relations for activators
#   1    2    3    4    5 NA's 
#   6    2    2   18    5    9
summary(as.factor(Repressors$bool))
# histogram of bool relations for repressors
#   1    2    3    4    5 NA's 
#  17    5   40  119   86   98 
load("../jul1.NAfiltered.tf.RData")
load("../jun27.GSE2180.gsplom.RData")
summary(as.factor(BOOL.TF.F[upper.tri(BOOL.TF.F)]))
#     1      2      3      4      5      6      7 
# 30984   3415  20991 106656  51782     85    272
sum(upper.tri(BOOL.TF.F)) # [1] 214185
summary(as.factor(BOOL.TF.F[upper.tri(BOOL.TF.F)]))/214185*100
##           1           2           3           4           5           6 
## 14.46599902  1.59441604  9.80040619 49.79620422 24.17629619  0.03968532 
##           7 
##  0.12699302 

sum(summary(as.factor(na.omit(Repressors$bool))))
#[1] 267
# fractions
#        1         2         3         4         5 
# 6.367041  1.872659 14.981273 44.569288 32.209738
sum(summary(as.factor(na.omit(Activators$bool[!self.act]))))
# 33
summary(as.factor(na.omit(Activators$bool[!self.act])))/33
#         1         2         3         4         5 
# 18.181818  6.060606  6.060606 54.545455 15.151515 
# some enrichment for MX in repressor, some un-enrichment for MX in activator


## weak? dcor?
summary(as.factor(Activators$weak[!self.act]))
## 2    3    5 NA's 
## 7   19    7    9
summary(as.factor(Repressors$weak))
##   1    2    3    4    5 NA's 
##   3   62  112   13   77   98 
WEAK.TF.F[upper.tri(WEAK.TF.F)]
##     1     2     3     4     5 
## 22620 95099 50589  9456 36421 
summary(as.factor(WEAK.TF.F[upper.tri(WEAK.TF.F)]))/sum(upper.tri(WEAK.TF.F))*100
#         1         2         3         4         5 
# 10.560964 44.400402 23.619301  4.414875 17.004459 
sum(summary(as.factor(na.omit(Repressors$weak))))
# [1] 267
#        1         2         3         4         5 
# 1.123596 23.220974 41.947566  4.868914 28.838951
# target detection necessary for row detection
# un-enriched for "and"

# Dcor
hist(DCOR.TF.F[upper.tri(DCOR.TF.F)])
hist(Repressors$dcor) # approx same dist
hist(Activators$dcor) # some skew to right


# ==================================================
# per-knockout map on GSPLOM

# GET GSPLOM
source("~/Dropbox/biostat/git_repos/dependency_glyph_splom/nec.net.plots.R")
source("~/Dropbox/biostat/git_repos/dependency_glyph_splom/lib.R")
# k=420
pdf("~/Desktop/GSE2180.gsplom.tf.filt.crits.jul15.k420.pdf", width=100, height=100)
R.GSE2180.F.TF <- splom(BOOL.TF.F, DCOR.TF.F, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.F.D)
z420 <- cutree(as.hclust(R.GSE2180.F.TF$Rhclust),k=420)
zz420 <- split(names(z420),z420)
syms <- rownames(DCOR.TF.F)[order.dendrogram(R.GSE2180.F.TF$Rhclust)]
R.z420.sig <- collapse.cls(BOOL.TF.F, z420, DCOR.TF.F, dcor.sig=0.36)
SUM.FLAWS <- R.z420.sig$MIX.SIGN + R.z420.sig$MIX.DIR
CRIT <- (log2(R.z420.sig$SIZE) <= SUM.FLAWS) & (SUM.FLAWS > 1)
clust.coords.420 <- clust.names.to.idx(zz420,syms)
plot.rects.coords(clust.coords.420)
map.plot.rects(clust.coords.420, CRIT&upper.tri(CRIT))
print(sum(CRIT&upper.tri(CRIT)))
dev.off()

save(R.GSE2180.F.TF, file="../jul30.k420.R.GSE2180.F.TF.RData")

# plot dendrogram
pdf("~/Desktop/GSE2180.gsplom.tf.filt.dendro.pdf", width=200, height=12)
plot(R.GSE2180.F.TF$Rhclust)
dev.off()



intersect(Repressors$entrez,Activators$entrez)
#integer(0)

gene <- "pal-1"
gsplom.order <- order.dendrogram(R.GSE2180.F.TF$Rhclust)

plot.reg.map <- function(gene, map.order, M) {
  gsplom.entrez <- colnames(M)[map.order]
  gsplom.syms <- rownames(M)[map.order]
  n <- length(gsplom.syms)
  repress.map <- gsplom.entrez %in% Repressors$entrez[Repressors$knockout==gene]
  activate.map <- gsplom.entrez %in% Activators$entrez[Activators$knockout==gene]
  pdf(paste0("~/Desktop/",gene,"regmap.pdf"), width=n/10, height=0.2)
  par(mar=c(0,0,0,0))
  image(1:n, 1:1, col=c("grey","red"), as.matrix(repress.map), axes=F)
  image(1:n, 1:1, col=c(rgb(0,0,0,0),"green"), as.matrix(activate.map), add=T, axes=F)
  axis(2, 1:n, labels=gsplom.syms, las = 2, line = -0.5, tick = 0, cex.axis = 0.7)
  dev.off()
}

for(g in unique(Tab.ko.edges$knockout)) {
  plot.reg.map(g, gsplom.order, DCOR.TF.F)
}

# TODO:
# pcc order for pal-1 reg targets
# z-score heatmap with reg targets

# map labeled matrix of activators, repressors


#gene.names <- unique(Tab.ko.edges$knockout)
vv <- sort(summary(as.factor(Activators$knockout)), decreasing=T)
gene.names <- names(vv)
act.cnts <- as.vector(vv)
# filter redundant tbx-9
act.cnts <- act.cnts[gene.names!="tbx-9"]
gene.names <- gene.names[gene.names!="tbx-9"]
m <- length(gene.names)
reg.map <- mat.or.vec(m,n)
for (i in 1:m) {
  gene <- gene.names[i]
  reg.map[i,gsplom.entrez %in% Repressors$entrez[Repressors$knockout==gene]] <- -1
  reg.map[i,gsplom.entrez %in% Activators$entrez[Activators$knockout==gene]] <- 1
  if (gene %in% Activators$genename[Activators$knockout==gene]) {
    reg.map[i,gsplom.syms==gene] <- 2
  } else {
    reg.map[i,gsplom.syms==gene] <- 3
  }
}

# merge tbx-8, tbx-9
# in GSPLOM order
pdf("~/Desktop/allgenes.regmap.pdf", width=3, height=30)
image(1:m, 1:n, col=c("red",rgb(0.8,0.8,0.8,1),"green","blue","black"), reg.map[,n:1], axes=F, breaks=c(-1,-0.5,0.5,1.5,2.5,3), xlab="", ylab="")
labs <- gene.names
labs[labs=="tbx-8"] <- "tbx-8*/9"
labs[labs %in% Activators$genename[Activators$knockout=="pal-1"]] <- paste0(labs[labs %in% Activators$genename[Activators$knockout=="pal-1"]],"*")
labs <- paste0(labs," (",act.cnts,")")
axis(1, 1:m, labels=labs, las = 2, line = -0.5, tick = 0, cex.axis = 0.7)
dev.off()


# get genes only in NA filtered matrix
qq <- match(rownames(BOOL.TF.F), rownames(M.TF))
M.TF.F <- M.TF[qq,]

# compute all-pairs PCC
library("gplots")
library("RColorBrewer")

CC <- cor(t(M.TF.F))
rownames(CC) <- rownames(DCOR.TF.F)
D.CC <- dist(CC)
CC.Rhclust <- as.dendrogram(hclust(D.CC, method="average"))
CC.Rowv <- rowMeans(CC, na.rm = TRUE)
CC.Rhclust <- reorder(CC.Rhclust, CC.Rowv)
heatmap.cols <- rev(colorRampPalette(brewer.pal(8,"RdBu"))(40))
heatmap.breaks <- seq(-1,1,0.05)

# TODO!
# Note: this is too unweildly. Manually draw the image using raster
pdf("~/Desktop/jul31.pcc.heatmap.pdf", width=80, height=80)
heatmap.2(CC, symm=TRUE, trace="none", Rowv=CC.Rhclust, Colv=CC.Rhclust, col=heatmap.cols, breaks=heatmap.breaks, key=F)
dev.off()
# make into AI figure


gene <- "pal-1"
pcc.heatmap.order <- order.dendrogram(CC.Rhclust)

plot.reg.map(gene, heatmap.order, CC)

### Knockout z-score heatmaps
load("../jul30.k420.R.GSE2180.F.TF.RData")
Tab.knockout.z <- read.table("../signed_zscore.csv", sep=",", header=T, stringsAsFactors=F)
gsplom.order <- order.dendrogram(R.GSE2180.F.TF$Rhclust)
gsplom.entrez <- colnames(BOOL.TF.F)[gsplom.order]
# pal-1
pal1.gsplom.z <- Tab.knockout.z$pal.1.z
pal1.gsplom.lg <- Tab.knockout.z$pal.1.lg
pal1.low.top <- (pal1.gsplom.z < -2) & (pal1.gsplom.lg > 1)
pal1.low <- (pal1.gsplom.z < -1) & (pal1.gsplom.lg > 0.3)
pal1.high.top <- (pal1.gsplom.z > 2) & (pal1.gsplom.lg < -1)
pal1.high <- (pal1.gsplom.z > 1) & (pal1.gsplom.lg < -0.3)

n <- length(pal1.low.top)
vv <- rep(0,n)
vv[pal1.low] <- -1
vv[pal1.low.top] <- -2
vv[pal1.high] <- 1
vv[pal1.high.top] <- 2
vv[Tab.knockout.z$pal.1.entrez=="175638"] <- 3

qq.pcc <- match(colnames(BOOL.TF.F)[order.dendrogram(CC.Rhclust)],Tab.knockout.z$pal.1.entrez)
qq.gsplom <- match(gsplom.entrez,Tab.knockout.z$pal.1.entrez)

pdf("~/Desktop/pcc.dendro.pdf", width=150, height=8)
plot(CC.Rhclust)
dev.off()


n<-655
pdf("~/Desktop/aug4.highlowmids.pdf", width=30, height=1)
par(mar=c(0,0,0,0))
image(1:n, 1:1, col=c("red","orange","white","#339900","#00ff00","#0000ff"), as.matrix(vv[qq.gsplom]), axes=F, xlab="", ylab="")
image(1:n, 1:1, col=c("red","orange","white","#339900","#00ff00","#0000ff"), as.matrix(vv[qq.pcc]), axes=F, xlab="", ylab="")
dev.off()
# save image
# compare with PCC order


library("gplots")
library("RColorBrewer")
heatmap.cols <- rev(colorRampPalette(brewer.pal(8,"RdBu"))(40))
heatmap.breaks <- c(-1000,seq(-4.75,4.75,0.25),1000)
pdf("~/Desktop/aug4.zscore.colheatmap.pdf", width=30, height=1)
par(mar=c(0,0,0,0))
image(1:length(pal1.gsplom.z), 1:1, col=heatmap.cols, breaks=heatmap.breaks, as.matrix(pal1.gsplom.z), xlab="", ylab="", axes=F)
dev.off()


# ------------------------------
# Aug 7th, 2013
# ------------------------------
# get top and marginal high and low KO results
Tab.knockout.z <- read.table("../signed_zscore.csv", sep=",", header=T, stringsAsFactors=F)
genes <- c("pal.1","lin.26","hnd.1","unc.120","nob.1","elt.3","nhr.25","hlh.1","tbx.8.9","elt.1","scrt.1","vab.7")
n <- length(Tab.knockout.z$pal.1.entrez)
gsplom.order <- order.dendrogram(R.GSE2180.F.TF$Rhclust)
gsplom.entrez <- colnames(BOOL.TF.F)[gsplom.order] # align entrez IDs to this mapping
gsplom.syms <- rownames(DCOR.TF.F)[gsplom.order]

to.order <- gsplom.order
entrez <- gsplom.entrez
syms <- gsplom.syms
high.lows <- list()
for (g in genes) {
  z <- Tab.knockout.z[[paste0(g,".z")]]
  lg <- Tab.knockout.z[[paste0(g,".lg")]]
  e <-  Tab.knockout.z[[paste0(g,".entrez")]]
  low.top <- (z < -2) & (lg > 1)
  low <- (z < -0.5) & (lg > 0.5)
  high.top <- (z > 2) & (lg < -1)
  high <- (z > 0.5) & (lg < -0.5)
  vv <- rep(0,n)
  vv[low] <- -1
  vv[low.top] <- -2
  vv[high] <- 1
  vv[high.top] <- 2
  vv <- vv[na.omit(match(entrez,e))]
  self <- sub("-",".",syms,fixed=T)==g
  if(length(vv[self])>0 && vv[self]==1) {
    vv[self] <- 4
  } else if (length(vv[self])>0 && vv[self]==2) {
    vv[self] <- 5
  } else {
    vv[self] <- 4
  }
  high.lows[[g]] <- vv
}

M <- mat.or.vec(length(entrez),length(genes))
for (i in 1:length(genes)) {
  M[,i] <- high.lows[[genes[i]]]
}
colnames(M) <- genes
rownames(M) <- syms

# high and self only
sort(apply(QQ>=2,2,sum),decreasing=T)
 ##  pal.1 tbx.8.9   elt.1   nob.1   hlh.1  scrt.1   hnd.1 unc.120   elt.3  nhr.25 
 ##      8       7       4       3       3       3       2       2       2       2 
 ## lin.26   vab.7 
 ##      1       1 



## pal-1 has most "strong" activators, but tbx-8/9 has most activators over-all (twice as many as pal-1)
## nhr-25, pal-1, and nob-1 are top three repressors by either threshold, tbx8/9 mid-ranked


# display in order of count of strong activators

strong.act.order <- c("pal.1","tbx.8.9","elt.1","hnd.1","scrt.1","hlh.1","nhr.25","nob.1","elt.3","lin.26","unc.120","vab.7")
M <- M[,match(strong.act.order,colnames(M))]
# append "gold" gene flag column
gold.genes <- c("pal-1", "tbx-8", "tbx-9", "elt-1", "cwn-1", "hnd-1", "scrt-1", "elt-3", "nhr-25", "nob-1", "vab-7", "unc-120", "hlh-1", "mab-21", "lin-26")
gold.flag <- as.numeric(syms %in% gold.genes)*3
gold.flag[gold.flag==0] <- -3
M <- cbind(gold.flag,M)
n<-length(entrez)
m<-length(genes)+1
pdf("~/Desktop/aug7.highlowmids.pdf", width=2, height=30)
par(mar=c(4,0,0,0))
image(1:m, 1:n, col=c("white","red","orange","#e5e5e5","#339900","#00ff00","#000000","#000099","#00ffff"), t(M[n:1,]), breaks=-3:6-0.5, axes=F, xlab="", ylab="")
axis(1, 1:m, labels=c("Gold Network",labs), las = 2, line = -0.5, tick = 0, cex.axis = 0.7)
dev.off()


pdf("activators.repressors.counts.pdf", width=12)
barplot(sort(apply(M>0&M!=4,2,sum),decreasing=T), main="Sum Activators")
barplot(sort(apply(M==2|M==5,2,sum),decreasing=T), main="Sum Strong Activators")
barplot(sort(apply(M<0,2,sum),decreasing=T), main="Sum Repressors")
barplot(sort(apply(M<=-2,2,sum),decreasing=T), main="Sum Strong Repressors")
dev.off()


pdf("negative.correlation.plots.pdf")
plot(M.TF[which(rownames(DCOR.TF)=="tbx-8"),], M.TF[which(rownames(DCOR.TF)=="hmg-3"),], main="sym neg")
plot(M.TF[which(rownames(DCOR.TF)=="tbx-8"),], M.TF[which(rownames(DCOR.TF)=="hmg-5"),])
plot(M.TF[which(rownames(DCOR.TF)=="hnd-1"),], M.TF[which(rownames(DCOR.TF)=="hmg-3"),], main="OR relation")
plot(M.TF[which(rownames(DCOR.TF)=="hnd-1"),], M.TF[which(rownames(DCOR.TF)=="hmg-5"),])
plot(M.TF[which(rownames(DCOR.TF)=="hmg-3"),], M.TF[which(rownames(DCOR.TF)=="hmg-5"),])
plot(M.TF[which(rownames(DCOR.TF)=="tbx-8"),], M.TF[which(rownames(DCOR.TF)=="hnd-1"),])
dev.off()


## generate network based on KO data
colnames(M) # nodes
rownames(M) %in% c(sub(".","-",colnames(M),fixed=T),"tbx-8",'tbx-9')
G <- M[which(rownames(M) %in% c(sub(".","-",colnames(M),fixed=T),"tbx-8",'tbx-9')),]
write.table(G, file="KO.adj.matrix.tab", sep="\t", col.names=NA)


## generate histograms for both strong and weaker relations
KO <- M
BOOL <- BOOL.TF.F[gsplom.order,gsplom.order]
WEAK <- WEAK.TF.F[gsplom.order,gsplom.order]
DCOR <- DCOR.TF.F[gsplom.order,gsplom.order]
rownames(BOOL) <- rownames(DCOR)
colnames(BOOL) <- rownames(DCOR)
rownames(WEAK) <- rownames(DCOR)
colnames(WEAK) <- rownames(DCOR)
colnames(DCOR) <- rownames(DCOR)

# exclude pal-1
ko.genes <- c("pal-1","hnd-1","unc-120","nob-1","elt-3","nhr-25","hlh-1","tbx-8","elt-1","scrt-1","tbx-9")
qq <- match(ko.genes,rownames(DCOR))
B <- BOOL[,qq]
ko.genes.col <- c("pal.1","hnd.1","unc.120","nob.1","elt.3","nhr.25","hlh.1","tbx.8.9","elt.1","scrt.1","tbx.8.9")
qq.col <- match(ko.genes.col,colnames(KO))
KO <- KO[,qq.col]
paste(colnames(B),colnames(KO)) # OK

all(rownames(B)==rownames(KO)) # true
summary(as.factor(B[KO %in% c(1,2)]))

sum(KO %in% c(1,2))
#[1] 680
sum(KO[,1:10] %in% c(1,2))
#[1] 571


summary(as.factor(B[KO %in% c(1,2)]))
#  1   2   3   4   5   6   7 
# 27   6  69 311 261   3   3
a <- as.list(summary(as.factor(B[KO %in% c(1,2)]))/680*100)
#         1          2          3          4          5          6          7 
# 3.9705882  0.8823529 10.1470588 45.7352941 38.3823529  0.4411765  0.4411765
n <- sum(upper.tri(BOOL.TF.F))*2 # 428370
c <- as.list(summary(as.factor(BOOL.TF.F[!diag(ncol(BOOL.TF.F))]))/n*100)
#           1           2           3           4           5           6 
# 12.13320261  1.59441604 12.13320261 49.79620422 24.17629619  0.03968532 
#           7 
#  0.12699302
sum(KO %in% c(-1,-2)) # 1383
b <- as.list(summary(as.factor(B[KO %in% c(-1,-2)])) /1383*100)
#          1           2           3           4           5           7 
# 8.89370933  3.10918294 12.65365148 51.12075199 24.15039769  0.07230658
b7 <- b[["7"]]
b[[6]] <- NULL
b[['6']] <- 0
b[['7']] <- b7


pdf("all.edges.hist.pdf", width=12)
barplot(rbind(as.numeric(a), as.numeric(b), as.numeric(c)), beside=T, xlab="Detection Level Class: Knockout (KO) versus Target (Y) lower thresh", ylab="Percent Total", names.arg=c("1: KO<=Y","2: KO<=>Y","3: KO=>Y","4: UNL","5: MX","6: NC","7: OR"), legend.text=c("Activating (680)", "Repressing (1383)", "All (428k)"))
dev.off()

# 

# tbx-8, hmg-3 time series plots

# KO plot for tbx-8, hmg-3
load("../celegans.apr8.expr.RData")
time <- phenoData(E.expr)$time[match(colnames(M.TF.F),colnames(E.expr))]
genotype <- phenoData(E.expr)$genotype[match(colnames(M.TF.F),colnames(E.expr))]
colors <- rep("black",length(genotype)) # NS: black
colors[genotype=="ms"] <- "red"         # ms: red

M.TF.F[rownames(DCOR.TF.F)=="tbx-8",]

library(celegansceentrezg.db)
ids <- rownames(E.expr)
GSE2180.syms <- mget(ids, celegansceentrezgSYMBOL, ifnotfound=NA)

pdf("~/Desktop/time.plots.pdf")
plot(time,M.TF.F[rownames(DCOR.TF.F)=="tbx-8",], col=colors)
plot(time,M.TF.F[rownames(DCOR.TF.F)=="elt-1",])
plot(time,M.TF.F[rownames(DCOR.TF.F)=="elt-3",])
plot(time,M.TF.F[rownames(DCOR.TF.F)=="nob-1",])
plot(time,M.TF.F[rownames(DCOR.TF.F)=="hlh-1",])
plot(time,M.TF.F[rownames(DCOR.TF.F)=="hnd-1",])
plot(time,M.TF.F[rownames(DCOR.TF.F)=="hmg-3",], col=colors)
plot(time,M.TF.F[rownames(DCOR.TF.F)=="pal-1",])
plot(time,M.TF.F[rownames(DCOR.TF.F)=="hmg-5",])
plot(time,M.TF.F[rownames(DCOR.TF.F)=="CELE_T13F2.2",], col=colors)
plot(time,M.TF.F[rownames(DCOR.TF.F)=="chd-1",], col=colors)

plot(time,M.TF.F[rownames(DCOR.TF.F)=="K09A11.1",], col=colors) # highly repressed
plot(time,M.TF.F[rownames(DCOR.TF.F)=="nhr-17",], col=colors)  # highly repressed


plot(time,M.TF.F[rownames(DCOR.TF.F)=="nhr-163",])
plot(time,M.TF.F[rownames(DCOR.TF.F)=="lin-32",], col=colors)


plot(time,M.TF.F[rownames(DCOR.TF.F)=="tbx-8",], col=colors) # all go up
plot(time,M.TF.F[rownames(DCOR.TF.F)=="lin-32",], col=colors) # only some go up
plot(M.TF.F[rownames(DCOR.TF.F)=="tbx-8",],M.TF.F[rownames(DCOR.TF.F)=="lin-32",], col=colors)
plot(M.TF.F[rownames(DCOR.TF.F)=="tbx-8",],M.TF.F[rownames(DCOR.TF.F)=="unc-130",], col=colors)
plot(M.TF.F[rownames(DCOR.TF.F)=="tbx-8",],M.TF.F[rownames(DCOR.TF.F)=="F21A10.2",], col=colors)
plot(M.TF.F[rownames(DCOR.TF.F)=="hmg-3",],M.TF.F[rownames(DCOR.TF.F)=="tbx-8",], col=colors)
plot(M.TF.F[rownames(DCOR.TF.F)=="hmg-3",],M.TF.F[rownames(DCOR.TF.F)=="hnd-1",], col=colors)
plot(M.TF.F[rownames(DCOR.TF.F)=="hmg-3",],M.TF.F[rownames(DCOR.TF.F)=="hlh-1",], col=colors)

plot(M.TF.F[rownames(DCOR.TF.F)=="tbx-8",],M.TF.F[rownames(DCOR.TF.F)=="nhr-17",], col=colors)
plot(M.TF.F[rownames(DCOR.TF.F)=="tbx-8",],M.TF.F[rownames(DCOR.TF.F)=="K09A11.1",], col=colors)
plot(M.TF.F[rownames(DCOR.TF.F)=="tbx-8",],M.TF.F[rownames(DCOR.TF.F)=="chd-1",], col=colors)


plot(M.TF.F[rownames(DCOR.TF.F)=="tbx-8",],M.TF.F[rownames(DCOR.TF.F)=="mex-5",], col=colors)
dev.off()

#mex plots
pdf("~/Desktop/mex.pie.pos.vs.tbx8.pdf")
plot(M.TF.F[rownames(DCOR.TF.F)=="tbx-8",],M.TF.F[rownames(DCOR.TF.F)=="mex-5",], col=colors, xlab="tbx-8", ylab="mex-5")
plot(M.TF.F[rownames(DCOR.TF.F)=="tbx-8",],M.TF.F[rownames(DCOR.TF.F)=="mex-6",], col=colors, xlab="tbx-8", ylab="mex-6")
plot(M.TF.F[rownames(DCOR.TF.F)=="tbx-8",],M.TF.F[rownames(DCOR.TF.F)=="mex-1",], col=colors, xlab="tbx-8", ylab="mex-1")
plot(M.TF.F[rownames(DCOR.TF.F)=="tbx-8",],M.TF.F[rownames(DCOR.TF.F)=="pos-1",], col=colors, xlab="tbx-8", ylab="pos-1")
plot(M.TF.F[rownames(DCOR.TF.F)=="tbx-8",],M.TF.F[rownames(DCOR.TF.F)=="pie-1",], col=colors, xlab="tbx-8", ylab="pie-1")
plot(time,M.TF.F[rownames(DCOR.TF.F)=="tbx-8",], col=colors, xlab="time", ylab="tbx-8")
plot(time,M.TF.F[rownames(DCOR.TF.F)=="mex-5",], col=colors, xlab="time", ylab="mex-5")
plot(time,M.TF.F[rownames(DCOR.TF.F)=="mex-6",], col=colors, xlab="time", ylab="mex-6")
plot(time,M.TF.F[rownames(DCOR.TF.F)=="mex-1",], col=colors, xlab="time", ylab="mex-1")
plot(time,M.TF.F[rownames(DCOR.TF.F)=="pos-1",], col=colors, xlab="time", ylab="pos-1")
plot(time,M.TF.F[rownames(DCOR.TF.F)=="pie-1",], col=colors, xlab="time", ylab="pie-1")
dev.off()


tbx8 <- M.TF.F[rownames(DCOR.TF.F)=="tbx-8",]
tbx9 <- M.TF.F[rownames(DCOR.TF.F)=="tbx-9",]
pdf("~/Desktop/tbx8.tbx9.pdf")
plot(time,tbx9, col=colors, xlab="time", ylab="tbx-9", ylim=c(0,3.5))
plot(time,tbx8, col=colors, xlab="time", ylab="tbx-8")
plot(tbx8,tbx9, col=colors, xlab="tbx-8", ylab="tbx-9", ylim=c(0,3.5))
plot(time,tbx8+tbx9, col=colors, xlab="time", ylab="tbx-8+tbx-9", ylim=c(0,4))
dev.off()



plot(M.TF.F[rownames(DCOR.TF.F)=="tbx-8",],M.TF.F[rownames(DCOR.TF.F)=="hmg-3",], col=colors, xlab="tbx-8", ylab="hmg-3")

#in mutuant: 1 is ms, 2 is N2

write.table(M.TF.F, file="../aug8.M.TF.F.tab", sep="\t", col.names=NA)
writeLines(as.character(genotype), "../aug8.genotype.tab")
writeLines(as.character(time), "../aug8.time.tab")

# add z-score, pcc, GSPLOM for time, mutant
# ---------------------------------------
# z-score
genotype <- as.factor(as.character(genotype))
split(M.TF.F[1,],genotype)

ZScore <- apply(M.TF.F, 1, function(v) t.test(v~genotype))
stats <- sapply(ZScore,function(v) v$statistic)
pvalues <- sapply(ZScore,function(v) v$p.value)
      

plot(time,M.TF.F[1,], col=colors, xlab="time", ylab=rownames(DCOR.TF.F)[1])


library("gplots")
library("RColorBrewer")o
heatmap.cols <- c('#0000ff',rev(colorRampPalette(brewer.pal(8,"RdBu"))(40)),'#ff0000') #42
heatmap.cols[21:22] <- rep('#ffffff',2)
heatmap.breaks <- c(-1000,seq(-5,5,0.25),1000)

pdf("~/Desktop/aug8.mutant.zstats.pdf", width=30, height=1)
par(mar=c(0,0,0,0))
image(1:length(stats), 1:1, col=heatmap.cols, breaks=heatmap.breaks, as.matrix(stats[gsplom.order]), xlab="", ylab="", axes=F)
dev.off()

rownames(DCOR.TF.F)[gsplom.order][1] # [1] "nhr-38"
stats[gsplom.order][1] # 3.379681 

pcc <- apply(M.TF.F,1,function(v)cor(v,time))

pdf("~/Desktop/aug8.time.pcc.pdf", width=30, height=1)
par(mar=c(0,0,0,0))
image(1:length(pcc), 1:1, col=heatmap.cols, breaks=c(-1,seq(-0.9,0.9,(1.8/40)),1), as.matrix(pcc[gsplom.order]), xlab="", ylab="", axes=F)
dev.off()

stats[rownames(DCOR.TF.F)=="tbx-8"] # -0.4543527 


plot(M.TF.F[rownames(DCOR.TF.F)=="tbx-8",],M.TF.F[rownames(DCOR.TF.F)=="tbx-2",], col=colors)

plot(M.TF.F[rownames(DCOR.TF.F)=="nhr-25",],M.TF.F[rownames(DCOR.TF.F)=="nhr-68",], col=colors)

plot(time,M.TF.F[rownames(DCOR.TF.F)=="rbp-10",], col=colors)

# no time pcc, high enrich for ms
plot(time,M.TF.F[rownames(DCOR.TF.F)=="rpb-10",], col=colors)

plot(time,M.TF.F[rownames(DCOR.TF.F)=="mbtr-1",], col=colors)
plot(time,M.TF.F[rownames(DCOR.TF.F)=="CELE_D1046.2",], col=colors)
plot(time,M.TF.F[rownames(DCOR.TF.F)=="mes-2",], col=colors)


plot(time,M.TF.F[rownames(DCOR.TF.F)=="nhr-222",], col=colors)
plot(time,M.TF.F[rownames(DCOR.TF.F)=="nhr-104",], col=colors)
plot(time,M.TF.F[rownames(DCOR.TF.F)=="npax-3",], col=colors)

plot(time,M.TF.F[rownames(DCOR.TF.F)=="tbx-33",], col=colors)
plot(time,M.TF.F[rownames(DCOR.TF.F)=="sdz-38",], col=colors)


plot(time,M.TF.F[rownames(DCOR.TF.F)=="nhr-163",], col=colors)

plot(time,M.TF.F[rownames(DCOR.TF.F)=="mis-2",], col=colors)


plot(time,M.TF.F[rownames(DCOR.TF.F)=="CELE_F23B12.7",], col=colors)
plot(time,M.TF.F[rownames(DCOR.TF.F)=="nhr-86",], col=colors)

K09A11.1

plot(time,M.TF.F[rownames(DCOR.TF.F)=="nhr-172",], col=colors)

plot(time,M.TF.F[rownames(DCOR.TF.F)=="F58G1.2",], col=colors)



pdf("~/Desktop/gold.clust.repress.outliers.pdf")
plot(time,M.TF.F[rownames(DCOR.TF.F)=="alr-1",], col=colors)
dev.off()

plot(time,M.TF.F[rownames(DCOR.TF.F)=="his-41",], col=colors)

# peak midway, enriched for pal-1 repress
plot(time,M.TF.F[rownames(DCOR.TF.F)=="F28C6.1",], col=colors)
plot(time,M.TF.F[rownames(DCOR.TF.F)=="egl-18",], col=colors)
plot(time,M.TF.F[rownames(DCOR.TF.F)=="CELE_F19F10.1",], col=colors)

T22C8.4


plot(time,M.TF.F[rownames(DCOR.TF.F)=="T22C8.4",], col=colors)
plot(time,M.TF.F[rownames(DCOR.TF.F)=="cnd-1",], col=colors)
plot(time,M.TF.F[rownames(DCOR.TF.F)=="lin-28",], col=colors)

# weird late activated even though mex-3 lower
plot(time,M.TF.F[rownames(DCOR.TF.F)=="nhr-256",], col=colors)
plot(time,M.TF.F[rownames(DCOR.TF.F)=="nhr-86",], col=colors)


# negative z score, pos cor with time, but not repressed
plot(time,M.TF.F[rownames(DCOR.TF.F)=="mep-1",], col=colors)
plot(time,M.TF.F[rownames(DCOR.TF.F)=="hda-1",], col=colors)
plot(time,M.TF.F[rownames(DCOR.TF.F)=="sex-1",], col=colors)
plot(time,M.TF.F[rownames(DCOR.TF.F)=="let-418",], col=colors)


# like mab-21, but N2 high rather than ms+N2 high
plot(time,M.TF.F[rownames(DCOR.TF.F)=="hlh-16",], col=colors)
plot(time,M.TF.F[rownames(DCOR.TF.F)=="ceh-10",], col=colors)



plot(time,M.TF.F[rownames(DCOR.TF.F)=="alr-1",], col=colors)
plot(time,M.TF.F[rownames(DCOR.TF.F)=="his-24",], col=colors)

plot(M.TF.F[rownames(DCOR.TF.F)=="alr-1",],M.TF.F[rownames(DCOR.TF.F)=="his-24",], col=colors)


plot(tbx8,M.TF.F[rownames(DCOR.TF.F)=="let-418",], col=colors)

plot(time,M.TF.F[rownames(DCOR.TF.F)=="ceh-32",], col=colors)
plot(time,M.TF.F[rownames(DCOR.TF.F)=="ceh-13",], col=colors)
plot(time,M.TF.F[rownames(DCOR.TF.F)=="elt-7",], col=colors)
plot(time,M.TF.F[rownames(DCOR.TF.F)=="ceh-36",], col=colors)
plot(time,M.TF.F[rownames(DCOR.TF.F)=="ceh-37",], col=colors)
plot(time,M.TF.F[rownames(DCOR.TF.F)=="ceh-43",], col=colors)
plot(time,M.TF.F[rownames(DCOR.TF.F)=="wrm-1",], col=colors)


#gold.genes
ems_clust <- c("med-2", "tbx-35", "ref-2", "alr-1")
cab_clust <- c("nhr-25", "nhr-171", "tbx-8", "tbx-9")
c_clust <- c("mab-5", "hlh-1", "pal-1", "cwn-1", "vab-7", "egl-5")
e_clust <- c("T28H10.3", "elt-2", "end-3", "elt-7", "dve-1", "end-1", "tps-2", "nhr-57", "nhr-79", "T23H4.2", "nhr-69", "nhr-68", "elt-2")

linage.color <- rep(0, length(gsplom.syms))
linage.color[gsplom.syms %in% gold.genes] <- 1
linage.color[gsplom.syms %in% c_clust] <- 1
linage.color[gsplom.syms %in% ems_clust] <- 2
linage.color[gsplom.syms %in% cab_clust] <- 3
linage.color[gsplom.syms %in% e_clust] <- 4

plot(time,M.TF.F[rownames(DCOR.TF.F)=="end-3",], col=colors)
plot(time,M.TF.F[rownames(DCOR.TF.F)=="nhr-69",], col=colors)

pdf("~/Desktop/aug10.linagemap.pdf", width=30, height=1)
par(mar=c(0,0,0,0))
image(1:length(gsplom.syms), 1:1, col=c("white","black","red","#666666","blue"), breaks=0:5-0.5, as.matrix(linage.color), xlab="", ylab="", axes=F)
dev.off()

# limited linage map
linage.color.mini <- rep(0, length(gsplom.syms))
linage.color.mini[gsplom.syms %in% gold.genes] <- 1
linage.color.mini[gsplom.syms %in% c_clust] <- 1
linage.color.mini[gsplom.syms %in% cab_clust] <- 1
linage.color.mini[gsplom.syms %in% e_clust] <- 2

pdf("~/Desktop/aug10.gold,c,abc_black--e_yellow.pdf", width=30, height=1)
par(mar=c(0,0,0,0))
image(1:length(gsplom.syms), 1:1, col=c("white","black","yellow"), breaks=0:3-0.5, as.matrix(linage.color.mini), xlab="", ylab="", axes=F)
dev.off()

plot(time,M.TF.F[rownames(DCOR.TF.F)=="end-1",], col=colors)


## ------------------------------
# Regulatory Inference
## ------------------------------
# As with all modeling, success in transcriptional regulatory inference from logical classes can depend on many known and unknown conditions, and to produce interpretable results can be challenging. To ground our inferences in biological relevance, in our methodology, we start from what we know, a set of well-understood genes that share a known relation. From these, we radiate outward to infer how these genes first relate to each other and then to other genes. In this example, the set of genes is chosen a priori by previous biological literature and the known relation is "positive regulation in the pal-1 dependent C linage." We call these pre-picked genes "landmarks" because they help us navigate an unknown landscape of relationships by orienting us with what is familiar.

# Because landmarks share a known relation, we assume that we can control for this relation by considering logical classes between landmarks only. That is, given that X and Y are both landmarks, if X is necessary for Y, we infer X regulates Y. [section on performance of this from BioVis paper]

# To infer regulatory relations between other transcription factors, we apply a strategy of "cluster and mark." First, we cluster genes using GSPLOM, then we identify clusters including landmarks. Because we assume that genes with similar expression patterns are related, we predict that genes clustered with landmarks share the landmark relation. Indeed, in support of this assumption, we observe that our pal-1  landmarks cluster together which highlighted in green on the GSPLOM. (Certainly other landmarks share the landmark relation!) By extension, we infer that the other genes clustered with our landmarks are also positively regulated in the pal-1 dependent transcription network.

# [possible sentence about how our method works but that alternative methods like differential expression of pal-1 linage-only over mixed linage misses important landmarks like tbx-8 (slightly negatively enriched in mex-3 mutants) and highlights a cluster that includes no landmarks.]

# Regulatory inference between and within GSPLOM clusters that do not include landmarks can be difficult because common relations in clusters without landmarks are unknown. However, with some knowledge about conditions represented in the data, we can make some regulatory inferences by using landmarks as references. For example, we know that the data modeled by the GSPLOM is a mixture of C linage only and mixed linage samples that include C linage. Also, from biological experience, we also know that transcription factors relevant to a lineage specification tend to start low and then increase in abundance over time as cells differentiate. Thus, we predict that pal-1 _negatively_ regulated genes are in non-landmark clusters, they are positively correlated with landmarks because they are related to linage specification, and they are not expressed in C linage only samples. That is, landmark gene high is asymmetrically necessary for possibly-negatively-regulated-gene X high. That is, in mixed linages, both the landmark and the prospective negative pal-1 target are highly expressed, but in C linage only, only the landmark is high.

#By figure [not included], we see that this is a reasonable inference. Either genes that meet these criteria are repressed by pal-1 or tbx-8/9 in KO experiements, or they are known to be specific to an alternate linage. Notably, we were able to make these inferences without knowing which samples belonged to which population or when during differentiation they were measured, only that meaningful variation in these dimensions were represented in the data. This is important as high quality homogenous, time series data are required by other modeling techiniques are difficult to acquire. In contrast, our method can be applied on more easily aquired unordered mixtures of samples, for example, large collections of biopsies of+ cancer mixed with different ratios of surrounding healthy tissue over many stages of progression.



Not negative correlation

pal-1 target necessary for other gene. Genes that are not expressed in C linage only

We also predict that nearly all genes with the landmark relation, pal-1 positively regulated genes relevant to the C linage, are clustered with landmarks. 

Thus, we can predict that pal-1 _negatively_ regulated genes are in another cluster

# Given clusters that include landmarks, we can make regulatory inferences in relation

We present one such strategy given that the data used to infer logical classes is known be to measure a mixture of pal-1 linage-only samples and mixed linage samples.



# Regulatory inference between two arbitrary genes can be difficult to model or interpret without 

# find clusters with known landmarks,
#   infer regulatory relations within these clusters / between landmarks

# inter-cluster, find relations that:
#   do not cluster with a known positive landmark
#   positive target high is necessary for putative repressed target high

# to form landmarked clusters, we merge together adjacent dendrogram branches
#  until such a merge would include an insignficant or a negative relation
