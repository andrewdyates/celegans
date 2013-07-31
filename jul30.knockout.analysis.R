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
integer(0)

gene <- "pal-1"

gsplom.order <- order.dendrogram(R.GSE2180.F.TF$Rhclust)
gsplom.entrez <- rownames(BOOL.TF.F)[gsplom.order]
gsplom.syms <- rownames(DCOR.TF.F)[gsplom.order]
n <- length(gsplom.names)
plot.reg.map <- function(gene) {
repress.map <- gsplom.entrez %in% Repressors$entrez[Repressors$knockout==gene]
activate.map <- gsplom.entrez %in% Activators$entrez[Activators$knockout==gene]
pdf(paste0("~/Desktop/",gene,"regmap.pdf"), width=n/10, height=0.2)
par(mar=c(0,0,0,0))
image(1:n, 1:1, col=c("grey","red"), as.matrix(repress.map), axes=F)
image(1:n, 1:1, col=c(rgb(0,0,0,0),"green"), as.matrix(activate.map), add=T, axes=F)
dev.off()
}

for(g in unique(Tab.ko.edges$knockout)) {
  plot.reg.map(g)
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
