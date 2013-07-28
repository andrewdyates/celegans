# local laptop paths
load("../jun27.GSE2180.gsplom.RData")
source("~/Dropbox/biostat/git_repos/dependency_glyph_splom/lib.R")

plot.rects <- function(zz, syms) {
  n <- length(syms)
  for(i in 1:length(zz)) {
    print(paste("cluster",i,"size",length(zz[[i]])))
    select.i <- which(syms %in% zz[[i]])
    x0 <- min(select.i)
    x1 <- max(select.i)
    y0 <- n-max(select.i)+1
    y1 <- n-min(select.i)+1
    rect(x0,y0,x1,y1, col=rgb(0,0,0,0.4))
  }
}

pdf("~/Desktop/GSE2180.gsplom.tf.jul1.k20.pdf", width=100, height=100)
R.GSE2180.TF <- splom(BOOL.TF, DCOR.TF, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.D)
n <- dim(BOOL.TF)[1]
z <- cutree(as.hclust(R.GSE2180.TF$Rhclust),k=20)
zz <- split(names(z),z)
syms <- rownames(DCOR.TF)[order.dendrogram(R.GSE2180.TF$Rhclust)]
plot.rects(zz)
dev.off()

pdf("~/Desktop/GSE2180.gsplom.tf.jul1.k30.pdf", width=100, height=100)
R.GSE2180.TF <- splom(BOOL.TF, DCOR.TF, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.D)
n <- dim(BOOL.TF)[1]
z <- cutree(as.hclust(R.GSE2180.TF$Rhclust),k=30)
zz <- split(names(z),z)
syms <- rownames(DCOR.TF)[order.dendrogram(R.GSE2180.TF$Rhclust)]
plot.rects(zz)
dev.off()

pdf("~/Desktop/GSE2180.gsplom.tf.jul1.k40.pdf", width=100, height=100)
R.GSE2180.TF <- splom(BOOL.TF, DCOR.TF, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.D)
n <- dim(BOOL.TF)[1]
z <- cutree(as.hclust(R.GSE2180.TF$Rhclust),k=40)
zz <- split(names(z),z)
syms <- rownames(DCOR.TF)[order.dendrogram(R.GSE2180.TF$Rhclust)]
plot.rects(zz)
dev.off()

pdf("~/Desktop/GSE2180.gsplom.tf.jul1.k50.pdf", width=100, height=100)
R.GSE2180.TF <- splom(BOOL.TF, DCOR.TF, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.D)
n <- dim(BOOL.TF)[1]
z <- cutree(as.hclust(R.GSE2180.TF$Rhclust),k=50)
zz <- split(names(z),z)
syms <- rownames(DCOR.TF)[order.dendrogram(R.GSE2180.TF$Rhclust)]
plot.rects(zz)
dev.off()

pdf("~/Desktop/GSE2180.gsplom.tf.jul1.k60.pdf", width=100, height=100)
R.GSE2180.TF <- splom(BOOL.TF, DCOR.TF, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.D)
n <- dim(BOOL.TF)[1]
z <- cutree(as.hclust(R.GSE2180.TF$Rhclust),k=60)
zz <- split(names(z),z)
syms <- rownames(DCOR.TF)[order.dendrogram(R.GSE2180.TF$Rhclust)]
plot.rects(zz)
dev.off()

pdf("~/Desktop/GSE2180.gsplom.tf.jul1.k70.pdf", width=100, height=100)
R.GSE2180.TF <- splom(BOOL.TF, DCOR.TF, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.D)
n <- dim(BOOL.TF)[1]
z <- cutree(as.hclust(R.GSE2180.TF$Rhclust),k=70)
zz <- split(names(z),z)
syms <- rownames(DCOR.TF)[order.dendrogram(R.GSE2180.TF$Rhclust)]
plot.rects(zz)
dev.off()

pdf("~/Desktop/GSE2180.gsplom.tf.jul1.k80.pdf", width=100, height=100)
R.GSE2180.TF <- splom(BOOL.TF, DCOR.TF, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.D)
n <- dim(BOOL.TF)[1]
z <- cutree(as.hclust(R.GSE2180.TF$Rhclust),k=80)
zz <- split(names(z),z)
syms <- rownames(DCOR.TF)[order.dendrogram(R.GSE2180.TF$Rhclust)]
plot.rects(zz)
dev.off()

pdf("~/Desktop/GSE2180.gsplom.tf.jul1.k10.pdf", width=100, height=100)
R.GSE2180.TF <- splom(BOOL.TF, DCOR.TF, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.D)
n <- dim(BOOL.TF)[1]
z <- cutree(as.hclust(R.GSE2180.TF$Rhclust),k=10)
zz <- split(names(z),z)
syms <- rownames(DCOR.TF)[order.dendrogram(R.GSE2180.TF$Rhclust)]
plot.rects(zz)
dev.off()

pdf("~/Desktop/GSE2180.gsplom.tf.dendro.pdf", width=180, height=12)
plot(R.GSE2180.TF$Rhclust, main="GSE2180 TF")
dev.off()

# Remove all self non-sym PC features
sum(diag(BOOL.TF)==0) # 13
sum(diag(BOOL.TF)==4) # 4
sum(diag(BOOL.TF)==2) # 742
# remove 17 features
pc.qq <- diag(BOOL.TF) == 2
BOOL.TF.F <- BOOL.TF[pc.qq,pc.qq]
DCOR.TF.F <- DCOR.TF[pc.qq,pc.qq]
WEAK.TF.F <- WEAK.TF[pc.qq,pc.qq]
# look for any 0 relations, remove offending features
sum(BOOL.TF==0)
sum(BOOL.TF.F==0)
sum.zeros <- apply(BOOL.TF.F==0,1,sum)

# strip out features with insufficient variances
while(sum(BOOL.TF.F==0)>0) {
  sum.zeros <- apply(BOOL.TF.F==0,1,sum)
  worst <- which.max(sum.zeros)
  qq <- 1:dim(BOOL.TF.F)[1]!=worst
  print(rownames(DCOR.TF.F)[!qq])
  BOOL.TF.F <- BOOL.TF.F[qq,qq]
  DCOR.TF.F <- DCOR.TF.F[qq,qq]
  WEAK.TF.F <- WEAK.TF.F[qq,qq]
  print(sum(BOOL.TF.F==0))
}
dim(BOOL.TF.F)
# 655

write.table(BOOL.TF.F, sep="\t", file="~/Desktop/jul2.GSE2180.bool.filt.tab")
BOOL.TF.F.D <- as.dist(as.matrix(read.table("~/Desktop/jul2.GSE2180.bool.filt.tab.booldist.tab", sep="\t", row.names=1, header=TRUE)))

pdf("~/Desktop/GSE2180.gsplom.tf.filt.jul2.pdf", width=100, height=100)
R.GSE2180.F.TF <- splom(BOOL.TF.F, DCOR.TF.F, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.F.D)
dev.off()


pdf("~/Desktop/GSE2180.gsplom.tf.filt.jul2.k20.pdf", width=100, height=100)
R.GSE2180.F.TF <- splom(BOOL.TF.F, DCOR.TF.F, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.F.D)
n <- dim(BOOL.TF.F)[1]
z <- cutree(as.hclust(R.GSE2180.F.TF$Rhclust),k=20)
zz <- split(names(z),z)
syms <- rownames(DCOR.TF.F)[order.dendrogram(R.GSE2180.F.TF$Rhclust)]
plot.rects(zz)
dev.off()


pdf("~/Desktop/GSE2180.gsplom.tf.filt.jul2.k30.pdf", width=100, height=100)
R.GSE2180.F.TF <- splom(BOOL.TF.F, DCOR.TF.F, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.F.D)
n <- dim(BOOL.TF.F)[1]
z <- cutree(as.hclust(R.GSE2180.F.TF$Rhclust),k=30)
zz <- split(names(z),z)
syms <- rownames(DCOR.TF.F)[order.dendrogram(R.GSE2180.F.TF$Rhclust)]
plot.rects(zz)
dev.off()


pdf("~/Desktop/GSE2180.gsplom.tf.filt.jul2.k40.avg.pdf", width=100, height=100)
R.GSE2180.F.TF.avg <- splom(BOOL.TF.F, DCOR.TF.F, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.F.D)
n <- dim(BOOL.TF.F)[1]
z <- cutree(as.hclust(R.GSE2180.F.TF.avg$Rhclust),k=40)
zz <- split(names(z),z)
syms <- rownames(DCOR.TF.F)[order.dendrogram(R.GSE2180.F.TF.avg$Rhclust)]
plot.rects(zz)
dev.off()

pdf("~/Desktop/GSE2180.gsplom.tf.filt.jul3.k40.complete.pdf", width=100, height=100)
R.GSE2180.F.TF.comp <- splom(BOOL.TF.F, DCOR.TF.F, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.F.D, clust.meth="complete")
n <- dim(BOOL.TF.F)[1]
z <- cutree(as.hclust(R.GSE2180.F.TF.comp$Rhclust),k=40)
zz <- split(names(z),z)
syms <- rownames(DCOR.TF.F)[order.dendrogram(R.GSE2180.F.TF.comp$Rhclust)]
plot.rects(zz)
dev.off()


# 1) find flaws, sum flaws, sign flaws
# 2) before and after splitting "loser" clusters with insufficient co-expression
#R.avg <- get.compression(BOOL.TF.F, as.hclust(R.GSE2180.F.TF.avg$Rhclust), DCOR.TF.F, min.dcor=0.36, max.k=300)
R.comp <- get.compression(BOOL.TF.F, as.hclust(R.GSE2180.F.TF.comp$Rhclust), DCOR.TF.F, min.dcor=0.36, max.k=300)


pdf("~/Desktop/GSE2180.gsplom.tf.filt.jul3.k40.pdf", width=100, height=100)
R.GSE2180.F.TF <- splom(BOOL.TF.F, DCOR.TF.F, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.F.D)
n <- dim(BOOL.TF.F)[1]
z <- cutree(as.hclust(R.GSE2180.F.TF$Rhclust),k=40)
zz <- split(names(z),z)
syms <- rownames(DCOR.TF.F)[order.dendrogram(R.GSE2180.F.TF$Rhclust)]
plot.rects(zz,syms)
dev.off()


# test at 300
H <- as.hclust(R.GSE2180.F.TF$Rhclust)
idx <- cutree(H,300)
C <- collapse.cls(BOOL.TF.F,idx,DCOR.TF.F)
S <- get.coh.M.score(C, 0.36)


pdf("~/Desktop/GSE2180.gsplom.tf.filt.jul3.k83.pdf", width=100, height=100)
R.GSE2180.F.TF <- splom(BOOL.TF.F, DCOR.TF.F, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.F.D)
n <- dim(BOOL.TF.F)[1]
z83 <- cutree(as.hclust(R.GSE2180.F.TF$Rhclust),k=83)
zz83 <- split(names(z83),z83)
syms <- rownames(DCOR.TF.F)[order.dendrogram(R.GSE2180.F.TF$Rhclust)]
plot.rects(zz83,syms)
dev.off()
R.z83 <- collapse.cls(BOOL.TF.F, z83, DCOR.TF.F)
# draw collapsed splom
pdf("~/Desktop/GSE2180.gsplom.tf.filt.jul3.k83.clusts.pdf", width=50, height=50)
R.K83 <- splom(R$CLS,R$DCOR,reorder=TRUE,MIN=0.36,sym=T)
dev.off()

# no dCor filtering
sum(R.z83$MIX.SIGN[upper.tri(R.z83$MIX.SIGN,diag=F)]>=1)
# 348
sum(R.z83$MIX.DIR[upper.tri(R.z83$MIX.SIGN,diag=F)]>=1)
# 228

# draw full splom with flawed edges highlighed
clust.names.to.idx <- function(zz,syms) {
  R <- list()
  n <- length(syms)
  for(i in 1:length(zz)) {
    R[[i]] <- list()
    select.i <- which(syms %in% zz[[i]])
    R[[i]]$x0 <- min(select.i)
    R[[i]]$x1 <- max(select.i)
    R[[i]]$y0 <- n-max(select.i)+1
    R[[i]]$y1 <- n-min(select.i)+1
  }
  R
}
plot.rects.coords <- function(coords) {
  for(c in coords) {
    rect(c$x0,c$y0,c$x1,c$y1, col=rgb(0,0,0,0.4))
  }
}
map.plot.rects <- function(coords,map) {
  ZZZ <- which(map, arr.ind=TRUE)
  for(i in 1:dim(ZZZ)[1]) {
    rect(coords[[ZZZ[i,1]]]$x0, coords[[ZZZ[i,2]]]$y0, coords[[ZZZ[i,1]]]$x1, coords[[ZZZ[i,2]]]$y1, col=rgb(0.5,0,0,0.3))
    rect(coords[[ZZZ[i,2]]]$x0, coords[[ZZZ[i,1]]]$y0, coords[[ZZZ[i,2]]]$x1, coords[[ZZZ[i,1]]]$y1, col=rgb(0.5,0,0,0.3))
  }
}


pdf("~/Desktop/GSE2180.gsplom.tf.filt.jul4.k83.pdf", width=100, height=100)
R.GSE2180.F.TF <- splom(BOOL.TF.F, DCOR.TF.F, useRaster=F, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.F.D)
z83 <- cutree(as.hclust(R.GSE2180.F.TF$Rhclust),k=83)
zz83 <- split(names(z83),z83)
syms <- rownames(DCOR.TF.F)[order.dendrogram(R.GSE2180.F.TF$Rhclust)]
clust.coords.83 <- clust.names.to.idx(zz83,syms)
plot.rects.coords(clust.coords.83)
Z.MIX <- R.z83$MIX.SIGN>0 & upper.tri(R.z83$MIX.SIGN, diag=FALSE)
map.plot.rects(clust.coords.83,Z.MIX)
dev.off()

# add significance filter to class matrix when computing mix signs
# only consider edges that are on average significant
z83 <- cutree(as.hclust(R.GSE2180.F.TF$Rhclust),k=83)
R.z83.sig <- collapse.cls(BOOL.TF.F, z83, DCOR.TF.F, dcor.sig=0.36)
sum(R.z83.sig$MIX.SIGN[upper.tri(R.z83.sig$MIX.SIGN,diag=F)]>=1)
# 44
sum(R.z83.sig$MIX.DIR[upper.tri(R.z83.sig$MIX.DIR,diag=F)]>=1)
# 192


pdf("~/Desktop/GSE2180.gsplom.tf.filt.jul5.k83.mixsign.clssig.pdf", width=100, height=100)
R.GSE2180.F.TF <- splom(BOOL.TF.F, DCOR.TF.F, useRaster=F, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.F.D)
z83 <- cutree(as.hclust(R.GSE2180.F.TF$Rhclust),k=83)
zz83 <- split(names(z83),z83)
syms <- rownames(DCOR.TF.F)[order.dendrogram(R.GSE2180.F.TF$Rhclust)]
clust.coords.83 <- clust.names.to.idx(zz83,syms)
plot.rects.coords(clust.coords.83)
Z.MIX <- R.z83.sig$MIX.SIGN>0 & upper.tri(R.z83.sig$MIX.SIGN, diag=FALSE)
map.plot.rects(clust.coords.83,Z.MIX)
dev.off()

# FIXED!
# there is a problem with one of the edges. Visually, in the GSPLOM,
# a contradiction seems to exist, yet the new algorithm does not flag it. why?

# edge is in nob-1(33), ceh-86(27) cluster edge
#> R.z83.sig$MIX.SIGN[33,27]
#[1] 1  (correct)
#> R.z83.sig$MIX.SIGN[27,33]
#[1] 0  (?)

sum(R.z83.sig$MIX.SIGN[upper.tri(R.z83.sig$MIX.SIGN,diag=F)]>=1)
# 44
sum(R.z83.sig$MIX.DIR[upper.tri(R.z83.sig$MIX.DIR,diag=F)]>=1)
# 214


pdf("~/Desktop/GSE2180.gsplom.tf.filt.jul5.k83.mixdir.clssig.pdf", width=100, height=100)
R.GSE2180.F.TF <- splom(BOOL.TF.F, DCOR.TF.F, useRaster=F, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.F.D)
Z.MIX <- R.z83.sig$MIX.DIR>0 & upper.tri(R.z83.sig$MIX.DIR, diag=FALSE)
plot.rects.coords(clust.coords.83)
map.plot.rects(clust.coords.83,Z.MIX)
dev.off()

# which features have the most flaws?


z350 <- cutree(as.hclust(R.GSE2180.F.TF$Rhclust),k=350)
R.z350.sig <- collapse.cls(BOOL.TF.F, z350, DCOR.TF.F, dcor.sig=0.36)

sum(R.z350.sig$MIX.SIGN[upper.tri(R.z350.sig$MIX.SIGN,diag=F)]>=1)
# 17
sum(R.z350.sig$MIX.DIR[upper.tri(R.z350.sig$MIX.DIR,diag=F)]>=1)
# 165


z500 <- cutree(as.hclust(R.GSE2180.F.TF$Rhclust),k=500)
R.z500.sig <- collapse.cls(BOOL.TF.F, z500, DCOR.TF.F, dcor.sig=0.36)
z400 <- cutree(as.hclust(R.GSE2180.F.TF$Rhclust),k=400)
R.z400.sig <- collapse.cls(BOOL.TF.F, z400, DCOR.TF.F, dcor.sig=0.36)
z450 <- cutree(as.hclust(R.GSE2180.F.TF$Rhclust),k=450)
R.z450.sig <- collapse.cls(BOOL.TF.F, z450, DCOR.TF.F, dcor.sig=0.36)


sum(R.z400.sig$MIX.SIGN[upper.tri(R.z400.sig$MIX.SIGN,diag=F)]>=1)
# 14
sum(R.z400.sig$MIX.DIR[upper.tri(R.z400.sig$MIX.DIR,diag=F)]>=1)
# 111

sum(R.z450.sig$MIX.SIGN[upper.tri(R.z450.sig$MIX.SIGN,diag=F)]>=1)
# 7
sum(R.z450.sig$MIX.DIR[upper.tri(R.z450.sig$MIX.DIR,diag=F)]>=1)
# 56

sum(R.z500.sig$MIX.SIGN[upper.tri(R.z500.sig$MIX.SIGN,diag=F)]>=1)
# 4
sum(R.z500.sig$MIX.DIR[upper.tri(R.z500.sig$MIX.DIR,diag=F)]>=1)
# 35


pdf("~/Desktop/GSE2180.gsplom.tf.filt.jul5.k83.mixdir.clssig.pdf", width=100, height=100)
R.GSE2180.F.TF <- splom(BOOL.TF.F, DCOR.TF.F, useRaster=F, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.F.D)
Z.MIX <- R.z350.sig$MIX.DIR>0 & upper.tri(R.z350.sig$MIX.DIR, diag=FALSE)
plot.rects.coords(clust.coords.83)
map.plot.rects(clust.coords.83,Z.MIX)
dev.off()


# TODO:
# compute coherence using dCor significance, mean coherence using only significant edges
# attempt to identify and break problem clusters
#  top clusters by flaws
#  any low dCor elements in cluster
#  any "highly incoherent" edges (20% flaw density)
save(BOOL.TF.F, BOOL.TF.F.D, DCOR.TF.F, WEAK.TF.F, file="../jul1.NAfiltered.tf.RData")
