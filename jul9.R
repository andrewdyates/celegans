# TODO:
# (+) compute coherence using dCor significance, mean coherence using only significant edges
# attempt to identify and break problem clusters
#  top clusters by flaws
#  any low dCor elements in cluster
#  any "highly incoherent" edges (log2 flaws)
# local laptop paths
# ---
# 1) split any clusters with any low-dcor relations
# 2) split any HI edges (or try to merge them in if on diag)

# TODO: confirm x-sign and x-dir and flaw counts
#   add "critically flawed" detector

load("../jun27.GSE2180.gsplom.RData")
load("../jul1.NAfiltered.tf.RData")
source("~/Dropbox/biostat/git_repos/dependency_glyph_splom/lib.R")
source("~/Dropbox/biostat/git_repos/dependency_glyph_splom/nec.net.plots.R")

#Z <- remove.na.features(BOOL.TF)

pdf("~/Desktop/GSE2180.gsplom.tf.filt.jul9.pdf", width=100, height=100)
R.GSE2180.F.TF <- splom(BOOL.TF.F, DCOR.TF.F, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.F.D)
dev.off()


# plot clusters and crit flaws
pdf("~/Desktop/GSE2180.gsplom.tf.filt.crits.jul15.pdf", width=100, height=100)
R.GSE2180.F.TF <- splom(BOOL.TF.F, DCOR.TF.F, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.F.D)
z83 <- cutree(as.hclust(R.GSE2180.F.TF$Rhclust),k=83)
zz83 <- split(names(z83),z83)
syms <- rownames(DCOR.TF.F)[order.dendrogram(R.GSE2180.F.TF$Rhclust)]
R.z83.sig <- collapse.cls(BOOL.TF.F, z83, DCOR.TF.F, dcor.sig=0.36)
SUM.FLAWS <- R.z83.sig$MIX.SIGN + R.z83.sig$MIX.DIR
CRIT <- (log2(R.z83.sig$SIZE) <= SUM.FLAWS) & (SUM.FLAWS > 1)
clust.coords.83 <- clust.names.to.idx(zz83,syms)
plot.rects.coords(clust.coords.83)
map.plot.rects(clust.coords.83, CRIT&upper.tri(CRIT))
print(sum(CRIT&upper.tri(CRIT)))
dev.off()




# k=200
# > [1] 29
# > [1] 0.9664277
pdf("~/Desktop/GSE2180.gsplom.tf.filt.crits.jul15.k200.pdf", width=100, height=100)
R.GSE2180.F.TF <- splom(BOOL.TF.F, DCOR.TF.F, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.F.D)
z200 <- cutree(as.hclust(R.GSE2180.F.TF$Rhclust),k=200)
zz200 <- split(names(z200),z200)
syms <- rownames(DCOR.TF.F)[order.dendrogram(R.GSE2180.F.TF$Rhclust)]
R.z200.sig <- collapse.cls(BOOL.TF.F, z200, DCOR.TF.F, dcor.sig=0.36)
SUM.FLAWS <- R.z200.sig$MIX.SIGN + R.z200.sig$MIX.DIR
CRIT <- (log2(R.z200.sig$SIZE) <= SUM.FLAWS) & (SUM.FLAWS > 1)
clust.coords.200 <- clust.names.to.idx(zz200,syms)
plot.rects.coords(clust.coords.200)
map.plot.rects(clust.coords.200, CRIT&upper.tri(CRIT))
print(sum(CRIT&upper.tri(CRIT)))
dev.off()
S.200 <- get.coh.M.score(R.z200.sig)
print(S.200$tri.wavg)

# k=300
# > [1] 11
# > [1] 0.9765304
pdf("~/Desktop/GSE2180.gsplom.tf.filt.crits.jul15.k300.pdf", width=100, height=100)
R.GSE2180.F.TF <- splom(BOOL.TF.F, DCOR.TF.F, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.F.D)
z300 <- cutree(as.hclust(R.GSE2180.F.TF$Rhclust),k=300)
zz300 <- split(names(z300),z300)
syms <- rownames(DCOR.TF.F)[order.dendrogram(R.GSE2180.F.TF$Rhclust)]
R.z300.sig <- collapse.cls(BOOL.TF.F, z300, DCOR.TF.F, dcor.sig=0.36)
SUM.FLAWS <- R.z300.sig$MIX.SIGN + R.z300.sig$MIX.DIR
CRIT <- (log2(R.z300.sig$SIZE) <= SUM.FLAWS) & (SUM.FLAWS > 1)
clust.coords.300 <- clust.names.to.idx(zz300,syms)
plot.rects.coords(clust.coords.300)
map.plot.rects(clust.coords.300, CRIT&upper.tri(CRIT))
print(sum(CRIT&upper.tri(CRIT)))
dev.off()
S.300 <- get.coh.M.score(R.z300.sig)
print(S.300$tri.wavg)

# coherence to 98%? then manually fix critical flaws?
# k=350
# > [1] 3
# > [1] 0.9810457
pdf("~/Desktop/GSE2180.gsplom.tf.filt.crits.jul15.k350.pdf", width=100, height=100)
R.GSE2180.F.TF <- splom(BOOL.TF.F, DCOR.TF.F, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.F.D)
z350 <- cutree(as.hclust(R.GSE2180.F.TF$Rhclust),k=350)
zz350 <- split(names(z350),z350)
syms <- rownames(DCOR.TF.F)[order.dendrogram(R.GSE2180.F.TF$Rhclust)]
R.z350.sig <- collapse.cls(BOOL.TF.F, z350, DCOR.TF.F, dcor.sig=0.36)
SUM.FLAWS <- R.z350.sig$MIX.SIGN + R.z350.sig$MIX.DIR
CRIT <- (log2(R.z350.sig$SIZE) <= SUM.FLAWS) & (SUM.FLAWS > 1)
clust.coords.350 <- clust.names.to.idx(zz350,syms)
plot.rects.coords(clust.coords.350)
map.plot.rects(clust.coords.350, CRIT&upper.tri(CRIT))
print(sum(CRIT&upper.tri(CRIT)))
dev.off()
S.350 <- get.coh.M.score(R.z350.sig)
print(S.350$tri.wavg)


# k=400
# > [1] 2
# > [1] 0.985512
pdf("~/Desktop/GSE2180.gsplom.tf.filt.crits.jul15.k400.pdf", width=100, height=100)
R.GSE2180.F.TF <- splom(BOOL.TF.F, DCOR.TF.F, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.F.D)
z400 <- cutree(as.hclust(R.GSE2180.F.TF$Rhclust),k=400)
zz400 <- split(names(z400),z400)
syms <- rownames(DCOR.TF.F)[order.dendrogram(R.GSE2180.F.TF$Rhclust)]
R.z400.sig <- collapse.cls(BOOL.TF.F, z400, DCOR.TF.F, dcor.sig=0.36)
SUM.FLAWS <- R.z400.sig$MIX.SIGN + R.z400.sig$MIX.DIR
CRIT <- (log2(R.z400.sig$SIZE) <= SUM.FLAWS) & (SUM.FLAWS > 1)
clust.coords.400 <- clust.names.to.idx(zz400,syms)
plot.rects.coords(clust.coords.400)
map.plot.rects(clust.coords.400, CRIT&upper.tri(CRIT))
print(sum(CRIT&upper.tri(CRIT)))
dev.off()
S.400 <- get.coh.M.score(R.z400.sig)
print(S.400$tri.wavg)

# k=450
# 0 crit flaws
# [1] 0.9896313
pdf("~/Desktop/GSE2180.gsplom.tf.filt.crits.jul15.k450.pdf", width=100, height=100)
R.GSE2180.F.TF <- splom(BOOL.TF.F, DCOR.TF.F, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.F.D)
z450 <- cutree(as.hclust(R.GSE2180.F.TF$Rhclust),k=450)
zz450 <- split(names(z450),z450)
syms <- rownames(DCOR.TF.F)[order.dendrogram(R.GSE2180.F.TF$Rhclust)]
R.z450.sig <- collapse.cls(BOOL.TF.F, z450, DCOR.TF.F, dcor.sig=0.36)
SUM.FLAWS <- R.z450.sig$MIX.SIGN + R.z450.sig$MIX.DIR
CRIT <- (log2(R.z450.sig$SIZE) <= SUM.FLAWS) & (SUM.FLAWS > 1)
clust.coords.450 <- clust.names.to.idx(zz450,syms)
plot.rects.coords(clust.coords.450)
map.plot.rects(clust.coords.450, CRIT&upper.tri(CRIT))
print(sum(CRIT&upper.tri(CRIT)))
dev.off()
S.450 <- get.coh.M.score(R.z450.sig)
print(S.450$tri.wavg)

# at what k is crit flaw count 0?
# k=425: 0
pdf("~/Desktop/GSE2180.gsplom.tf.filt.crits.jul15.k425.pdf", width=100, height=100)
R.GSE2180.F.TF <- splom(BOOL.TF.F, DCOR.TF.F, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.F.D)
z425 <- cutree(as.hclust(R.GSE2180.F.TF$Rhclust),k=425)
zz425 <- split(names(z425),z425)
syms <- rownames(DCOR.TF.F)[order.dendrogram(R.GSE2180.F.TF$Rhclust)]
R.z425.sig <- collapse.cls(BOOL.TF.F, z425, DCOR.TF.F, dcor.sig=0.36)
SUM.FLAWS <- R.z425.sig$MIX.SIGN + R.z425.sig$MIX.DIR
CRIT <- (log2(R.z425.sig$SIZE) <= SUM.FLAWS) & (SUM.FLAWS > 1)
clust.coords.425 <- clust.names.to.idx(zz425,syms)
plot.rects.coords(clust.coords.425)
map.plot.rects(clust.coords.425, CRIT&upper.tri(CRIT))
print(sum(CRIT&upper.tri(CRIT)))
dev.off()
S.425 <- get.coh.M.score(R.z425.sig)
print(S.425$tri.wavg)


# k=412: 2
pdf("~/Desktop/GSE2180.gsplom.tf.filt.crits.jul15.k412.pdf", width=100, height=100)
R.GSE2180.F.TF <- splom(BOOL.TF.F, DCOR.TF.F, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.F.D)
z412 <- cutree(as.hclust(R.GSE2180.F.TF$Rhclust),k=412)
zz412 <- split(names(z412),z412)
syms <- rownames(DCOR.TF.F)[order.dendrogram(R.GSE2180.F.TF$Rhclust)]
R.z412.sig <- collapse.cls(BOOL.TF.F, z412, DCOR.TF.F, dcor.sig=0.36)
SUM.FLAWS <- R.z412.sig$MIX.SIGN + R.z412.sig$MIX.DIR
CRIT <- (log2(R.z412.sig$SIZE) <= SUM.FLAWS) & (SUM.FLAWS > 1)
clust.coords.412 <- clust.names.to.idx(zz412,syms)
plot.rects.coords(clust.coords.412)
map.plot.rects(clust.coords.412, CRIT&upper.tri(CRIT))
print(sum(CRIT&upper.tri(CRIT)))
dev.off()
S.412 <- get.coh.M.score(R.z412.sig)
print(S.412$tri.wavg)

# k=419: 1
pdf("~/Desktop/GSE2180.gsplom.tf.filt.crits.jul15.k419.pdf", width=100, height=100)
R.GSE2180.F.TF <- splom(BOOL.TF.F, DCOR.TF.F, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.F.D)
z419 <- cutree(as.hclust(R.GSE2180.F.TF$Rhclust),k=419)
zz419 <- split(names(z419),z419)
syms <- rownames(DCOR.TF.F)[order.dendrogram(R.GSE2180.F.TF$Rhclust)]
R.z419.sig <- collapse.cls(BOOL.TF.F, z419, DCOR.TF.F, dcor.sig=0.36)
SUM.FLAWS <- R.z419.sig$MIX.SIGN + R.z419.sig$MIX.DIR
CRIT <- (log2(R.z419.sig$SIZE) <= SUM.FLAWS) & (SUM.FLAWS > 1)
clust.coords.419 <- clust.names.to.idx(zz419,syms)
plot.rects.coords(clust.coords.419)
map.plot.rects(clust.coords.419, CRIT&upper.tri(CRIT))
print(sum(CRIT&upper.tri(CRIT)))
dev.off()
S.419 <- get.coh.M.score(R.z419.sig)
print(S.419$tri.wavg)

# k=422: 0
pdf("~/Desktop/GSE2180.gsplom.tf.filt.crits.jul15.k422.pdf", width=100, height=100)
R.GSE2180.F.TF <- splom(BOOL.TF.F, DCOR.TF.F, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.F.D)
z422 <- cutree(as.hclust(R.GSE2180.F.TF$Rhclust),k=422)
zz422 <- split(names(z422),z422)
syms <- rownames(DCOR.TF.F)[order.dendrogram(R.GSE2180.F.TF$Rhclust)]
R.z422.sig <- collapse.cls(BOOL.TF.F, z422, DCOR.TF.F, dcor.sig=0.36)
SUM.FLAWS <- R.z422.sig$MIX.SIGN + R.z422.sig$MIX.DIR
CRIT <- (log2(R.z422.sig$SIZE) <= SUM.FLAWS) & (SUM.FLAWS > 1)
clust.coords.422 <- clust.names.to.idx(zz422,syms)
plot.rects.coords(clust.coords.422)
map.plot.rects(clust.coords.422, CRIT&upper.tri(CRIT))
print(sum(CRIT&upper.tri(CRIT)))
dev.off()
S.422 <- get.coh.M.score(R.z422.sig)
print(S.422$tri.wavg)

# k=420
#> [1] 0.9873064
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
S.420 <- get.coh.M.score(R.z420.sig)
print(S.420$tri.wavg)

# OK! 420 it is.
