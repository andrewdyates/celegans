library(Biobase)
load("../trans.genes.apr17.gsplom.RData")
load("../apr16.genelevel.depM.RData")
load("../apr16.genelevel.exprs.RData")
source("~/Dropbox/biostat/git_repos/dependency_glyph_splom/lib.R")

make.cls.dist.M <- function() {
  cls.enum.names <- c("hih","pc","lil","unl","hil","nc","lih","na")
  r1 <- c(0,1,2,1,2,3,2,2)
  r2 <- c(1,0,1,2,3,4,3,2)
  r3 <- c(2,1,0,1,2,3,2,2)
  r4 <- c(1,2,1,0,1,2,1,2)
  r5 <- c(2,3,2,1,0,1,2,2)
  r6 <- c(3,4,3,2,1,0,1,2)
  r7 <- c(2,3,2,1,2,1,0,2)
  r8 <- c(2,2,2,2,2,2,2,0)
  G.DIST.TABLE <- rbind(r1,r2,r3,r4,r5,r6,r7,r8)
  rownames(G.DIST.TABLE) <- cls.enum.names
  colnames(G.DIST.TABLE) <- cls.enum.names
  G.DIST.TABLE
}
G.DIST <- make.cls.dist.M()

# coherences is max 1-(1/4) mean hamming dist to single class
# BOOL_ENUM = {0:'NA', 1:'XiY', 2:'PC', 3:'YiX', 4:'UNL', 5:'MX', 6:'NC', 7:'OR'}
# NOTE: NA:8 in R CLS dist matrix
# Count occurences of classes. Indexed from 1. NA:0 is index 8.
get.cls.counts <- function(CLS, sym=T) {
  if (sym) {
    stopifnot(dim(CLS)[1] == dim(CLS)[2])
    C <- CLS[upper.tri(CLS)]
  } else {
    C <- c(CLS)
  }
  # indexed from 1. NA is index 8
  if(any(CLS==0)) stop("Set NA CLS==0 to 8.")
  cls.counts <- sapply(c(1:8), function(i) sum(C==i))
  # account for arbitrary symmetry across diagonal for XiY and YiX classes
  if (sym) {
    cls.counts[1] = cls.counts[1] + cls.counts[3]
    cls.counts[3] = 0
  }
  cls.counts
}
# Calculate coherence for a cls given a count of other classes.
#   cls.counts enumerated from 1. NA is index 8 (not index 0)
get.coh <- function(cls.counts, cls) {
  stopifnot(length(cls.counts) == 8)
  stopifnot(cls >= 1 && cls <= 8)
  n <- sum(cls.counts)
  sigma <- sum(cls.counts * G.DIST[cls,])
  e <- 1-(sigma/n/4)
  (e-0.5)*2 # center at 0, range -1 to 1
}
# Calculate coherence coherence vector for all classes
#   cls.counts enumerated from 1. NA is index 8 (not index 0)
get.coh.vec <- function(CLS, sym=T) {
  if (length(CLS) == 1) sym<-F
  cls.counts <- get.cls.counts(CLS, sym)
  sapply(1:8, function(i) get.coh(cls.counts, i))
}

# Return most coherent class enumeration given coherence vector.
#   coh.v enumerated from 1. NA is index 8 (not index 0)
choose.coh.cls <- function(coh.v, min.coh=0.74) {
  max.coh <- max(coh.v)
  max.coh.i <- coh.v == max.coh
  #print(coh.v)
  #print(max.coh.i)
  # coherence below threshold. Return
  R <- list(); R$coh <- max.coh
  b <- max.coh < min.coh
  if (is.na(b)) {
    print(coh.v)
    stop("BAD")
  }
  if (b) {
    R$cls <- 4
  } else {
    if (sum(max.coh.i)==1) {
      R$cls <- which(max.coh.i)
    # Handle ties.
    # 1. ignore NA class
    } else if (sum(max.coh.i[1:7])==1) {
      R$cls <- which(max.coh.i[1:7])
    # 1.5 ignore NA and UNL class
    } else if (sum(max.coh.i[c(1,2,3,5,6,7)]) == 1) {
      max.coh.i[4] <- FALSE
      R$cls <- which(max.coh.i[1:7])
    # 2. if ties for opposite signs, return UNL
    } else if (any(max.coh.i[c(1,2,3)]) & any(max.coh.i[c(5,6,7)])) {
      R$cls <- 4
    # 3. if all of a sign or both asyms of a sign, return sym of that sign
    } else if (all(max.coh.i[c(1,2,3)]) | all(max.coh.i[c(1,3)])) {
      R$cls <- 2
    } else if (all(max.coh.i[c(5,6,7)]) | all(max.coh.i[c(5,7)])) {
      R$cls <- 6
    # 4. if an asym and sym of the same sign, return asym
    } else if (all(max.coh.i[c(1,2)]) | all(max.coh.i[c(2,3)])) {
      R$cls <- 2
    } else if (all(max.coh.i[c(5,6)]) | all(max.coh.i[c(6,7)])) {
      R$cls <- 6
    }
  }
  R
}

get.mean.dcor <- function(DCOR, sym=F) {
  if (sym) {
    D <- DCOR[upper.tri(DCOR)]
  } else {
    D <- DCOR
  }
  u<-mean(D)
  if(is.na(u))
    u <- 1
  u
}

# Return collapsed class, coherence, and overall coherence by indices.
collapse.cls <- function(CLS, idx, DCOR=NULL) {
  CLS[CLS==0] <- 8
  idx <- as.factor(idx)
  n <- length(levels(idx))
  R <- list()
  R$idx <- idx
  R$SIZE <- matrix(0, nrow=n, ncol=n)
  R$CLS <- matrix(0, nrow=n, ncol=n)
  rownames(R$CLS) <- 1:n
  colnames(R$CLS) <- 1:n
  R$COH <- matrix(0, nrow=n, ncol=n)
  R$MIX.SIGN <- matrix(0.0, nrow=n, ncol=n)
  R$MIX.DIR <- matrix(0.0, nrow=n, ncol=n)
  R$members <- split(rownames(CLS), idx)
  if (!is.null(DCOR)) {
    R$DCOR <- matrix(0, nrow=n, ncol=n)
    rownames(R$DCOR) <- 1:n
    colnames(R$DCOR) <- 1:n
  } else {
    R$DCOR <- NULL
  }
  # compute overall coherence
  for (i in 1:length(levels(idx))) {
    for (j in 1:length(levels(idx))) {
      iv <- idx == levels(idx)[i]
      jv <- idx == levels(idx)[j]
      C <- CLS[iv,jv]
      if (!is.null(DCOR))
        R$DCOR[i,j] <- get.mean.dcor(DCOR[iv,jv], sym=i==j)
      r <- choose.coh.cls(get.coh.vec(C, sym=i==j))
      if (is.null(r$cls)) {
        print(C)
        print(r)
        stop("Null class.")
      }
      R$CLS[i,j] <- r$cls
      R$COH[i,j] <- r$coh
      R$MIX.SIGN[i,j] <- has.mix.sign(C)
      R$MIX.DIR[i,j]  <- has.mix.direction(C)
      if (i == j) {
        stopifnot(sum(iv)==sum(jv))
        R$SIZE[i,j] <- sum(iv)*(sum(iv)-1)/2
      } else {
        R$SIZE[i,j] <- sum(iv)*sum(jv)
      }
    }
  }
  R
}

# count number of cross edges
has.mix.sign <- function(CLS)
  any(c(1,2,3) %in% CLS) & any(c(5,6,7) %in% CLS)

has.mix.direction <- function(CLS)
  all(c(1,3) %in% CLS) | all(c(5,7) %in% CLS)

get.wavg.score <- function(coh, size) {
  nn <- sum(size)
  p <- size/nn
  sum(p*coh)
}

get.coh.M.score <- function(COLLAPSED, min.dcor=0) {
  R <- list()
  all.tri <- upper.tri(COLLAPSED$COH,diag=F)
  all.diag <- upper.tri(COLLAPSED$COH,diag=T) & !all.tri
  if (is.null(COLLAPSED$DCOR)) {
    tri <- all.tri
    diag <- all.diag
    R$dcor.clust <- NULL
  } else {
    d <- COLLAPSED$DCOR >= min.dcor
    tri <- d & all.tri
    diag <- d & all.diag
    R$min.dcor <- min.dcor
    R$num.dcor.clust <- sum(d)
  }
  R$all.wavg <- get.wavg.score(COLLAPSED$COH, COLLAPSED$SIZE)
  R$tri.wavg <- get.wavg.score(COLLAPSED$COH[tri], COLLAPSED$SIZE[tri])
  R$all.tri.wavg <- get.wavg.score(COLLAPSED$COH[all.tri], COLLAPSED$SIZE[all.tri])
  R$diag.wavg <- get.wavg.score(COLLAPSED$COH[diag], COLLAPSED$SIZE[diag])
  R$all.diag.wavg <- get.wavg.score(COLLAPSED$COH[all.diag], COLLAPSED$SIZE[all.diag])
  R$edge.xsign <- sum(COLLAPSED$MIX.SIGN[tri])
  R$edge.xdir <- sum(COLLAPSED$MIX.DIR[tri])
  R$edge.sum.flaws <- sum(COLLAPSED$MIX.SIGN[tri]) + sum(COLLAPSED$MIX.DIR[tri])
  z <- COLLAPSED$MIX.SIGN[tri] | COLLAPSED$MIX.DIR[tri]
  R$edge.flaws <- sum(z)
  R$edge.all.n <- sum(all.tri)
  R$edge.n <- sum(tri)
  R$clust.all.n <- sum(all.diag)
  R$clust.n <- sum(diag)
  R
}

# ignore contribution of edges below dCor threshold
get.compression <- function(CLS, H, DCOR=NULL, min.dcor=0) {
  #n <- round(dim(CLS)[1]/2)+2
  n <- dim(CLS)[1]
  R <- list()
  #R$clust.wavg <- rep(0,n)
  #R$edge.wavg <- rep(0,n)
  #R$edge.flaws <- rep(0,n)
  for(k in 1:n) {
    print(k)
    idx <- cutree(H,k)
    C <- collapse.cls(CLS,idx,DCOR)
    S <- get.coh.M.score(C, min.dcor)
    R$all.wavg[k] <- S$all.wavg
    R$clust.wavg[k] <- S$diag.wavg
    R$edge.wavg[k] <- S$tri.wavg
    
    R$all.clust.wavg[k] <- S$all.diag.wavg
    R$all.edge.wavg[k] <- S$all.tri.wavg

    R$edge.flaws[k] <- S$edge.flaws
    R$edge.sum.flaws[k] <- S$edge.sum.flaws
    R$edge.xsign[k] <- S$edge.xsign
    R$edge.xdir[k] <- S$edge.xdir
    R$edge.n[k] <- S$edge.n
    R$edge.all.n[k] <- S$edge.all.n
    R$clust.n[k] <- S$clust.n
    R$clust.all.n[k] <- S$clust.all.n
    if (!is.null(S$dcor.clust))
      R$dcor.clust[k] <- S$dcor.clust
  }
  R
}

## pdf("/dev/null")
## GOLD.R <- splom(D.expr.gold$CLS, D.expr.gold$DCOR); dev.off();
## H <- as.hclust(GOLD.R$Rhclust)
## CLS <- D.expr.gold$CLS
## DCOR <- D.expr.gold$DCOR

## idx <- cutree(H,k=5)
## R <- collapse.cls(CLS, idx, DCOR)
## get.coh.M.score(R, min.dcor=0.3)

## Z <- get.compression(CLS,H,DCOR,min.dcor=0.32)

## H <- as.hclust(TRANS.R$Rhclust)
## CLS <- D.expr.trans$CLS
## DCOR <- D.expr.trans$DCOR
## # RUN THIS
## ZZ <- get.compression(CLS, H, DCOR, min.dcor=0.32)


## plot(ZZ$edge.n / ZZ$edge.all.n)
## plot(ZZ$clust.wavg)

## plot(ZZ$edge.wavg, ZZ$all.tri.wavg)
## plot(ZZ$all.tri.wavg)

## plot(ZZ$edge.wavg[5:191], ZZ$all.tri.wavg[5:191])
## plot(ZZ$edge.flaws)
## plot(ZZ$edge.flaws / ZZ$edge.all.n)
## plot(ZZ$edge.n / ZZ$edge.all.n)
## plot(ZZ$edge.wavg)
## plot(ZZ$edge.wavg[5:191])

## ZZ$edge.flaws / ZZ$edge.all.n <0.05
## save(ZZ, file="quant.all.trans.k.RData")

# --------------------------------------------------
# k == 70
# edge weight > 90% for extant edges
# under 5% edge flaw for extant edges

stopifnot(all(order.dendrogram(TRANS.R$Rhclust) == order.dendrogram(TRANS.R$Chclust)))
H <- as.hclust(TRANS.R$Rhclust)
CLS <- D.expr.trans$CLS
DCOR <- D.expr.trans$DCOR
idx <- cutree(H,k=70)
R <- collapse.cls(CLS, idx, DCOR)

S <- get.coh.M.score(R, min.dcor=0.32)

#hist(sqrt(R$SIZE))

# GSPLOM BUG?
pdf("all.trans.gsplom.k70.d32.apr27.pdf", width=20, height=20)
G <- splom(R$CLS, R$DCOR, grid=F, MIN=0.32, MAX=1, sym=T)
dev.off()

pdf("all.trans.dendro.k70.d32.apr27.pdf", width=20, height=8)
plot(G$Rhclust)
dev.off()


## Attempt to compress compressed matrix
H <- as.hclust(G$Rhclust)
CLS <- R$CLS
DCOR <- R$DCOR
ZZZ <- get.compression(CLS,H,DCOR,min.dcor=0.32)


plot(ZZZ$edge.wavg[4:70])
plot(ZZZ$edge.flaws / ZZZ$edge.all.n)
plot(ZZZ$edge.n / ZZZ$edge.all.n)

th <- 0.05 / (dim(D.expr.trans$CLS)[1] / dim(R$CLS)[1])

# compressed flaws
plot(ZZZ$edge.flaws / ZZZ$edge.all.n)
flaw.ratio <- ZZZ$edge.flaws / ZZZ$edge.all.n
which(flaw.ratio <= th)[1] # 43

plot(ZZZ$edge.wavg)
which(ZZZ$edge.wavg >= 0.9)[1] # 44

## Compress k70 to k44
## ------------------------------
H <- as.hclust(G$Rhclust)
CLS <- R$CLS
DCOR <- R$DCOR
idx <- cutree(H,k=44)
RR <- collapse.cls(CLS, idx, DCOR)
SS <- get.coh.M.score(RR, min.dcor=0.32)


pdf("all.trans.gsplom.k70.k44.d32.apr27.pdf", width=20, height=20)
GG <- splom(RR$CLS, RR$DCOR, grid=F, MIN=0.32, MAX=1, sym=T)
dev.off()

pdf("all.trans.dendro.k70.k44.d32.apr27.pdf", width=20, height=8)
plot(GG$Rhclust)
dev.off()

## Can we compress k44?
H <- as.hclust(GG$Rhclust)
CLS <- RR$CLS
DCOR <- RR$DCOR
ZZZZ <- get.compression(CLS,H,DCOR,min.dcor=0.32)
th2 <- 0.05 / (dim(D.expr.trans$CLS)[1] / dim(RR$CLS)[1])

flaw.ratio2 <- ZZZZ$edge.flaws / ZZZZ$edge.all.n
plot(flaw.ratio2)
which(flaw.ratio2 <= th2) # 35
plot(ZZZZ$edge.wavg)
which(ZZZZ$edge.wavg >= 0.9)[1] # 31
## compressed to k=35!
idx <- cutree(H,k=35)
RRR <- collapse.cls(CLS, idx, DCOR)
SSS <- get.coh.M.score(RRR, min.dcor=0.32)

pdf("all.trans.gsplom.k70.k44.k35.d32.apr27.pdf", width=20, height=20)
GGG <- splom(RRR$CLS, RRR$DCOR, grid=F, MIN=0.32, MAX=1, sym=T)
dev.off()

pdf("all.trans.dendro.k70.k44.k35.d32.apr27.pdf", width=20, height=8)
plot(GGG$Rhclust)
dev.off()

## --------------------------------------------
## Can we compress k35?
H <- as.hclust(GGG$Rhclust)
CLS <- RRR$CLS
DCOR <- RRR$DCOR
ZZZZZ <- get.compression(CLS,H,DCOR,min.dcor=0.32)
th3 <- 0.05 / (dim(D.expr.trans$CLS)[1] / dim(RRR$CLS)[1])
flaw.ratio3 <- ZZZZZ$edge.flaws / ZZZZZ$edge.all.n
which(flaw.ratio3 <= th3)[1] # 31
which(ZZZZZ$edge.wavg >= 0.9)[1] # 27
## compressed k=35 to k=31!
idx <- cutree(H,k=31)
RRRR <- collapse.cls(CLS, idx, DCOR)
SSSS <- get.coh.M.score(RRRR, min.dcor=0.32)

pdf("all.trans.gsplom.k70.k44.k35.k31.d32.apr27.pdf", width=20, height=20)
GGGG <- splom(RRRR$CLS, RRRR$DCOR, grid=F, MIN=0.32, MAX=1, sym=T)
dev.off()

pdf("all.trans.dendro.k70.k44.k35.k31.d32.apr27.pdf", width=20, height=8)
plot(GGGG$Rhclust)
dev.off()

# ==================================================
## compressed k35 to k31. Can we compress k31?
H <- as.hclust(GGGG$Rhclust)
CLS <- RRRR$CLS
DCOR <- RRRR$DCOR
ZZZZZZ <- get.compression(CLS,H,DCOR,min.dcor=0.32)
th4 <- 0.05 / (dim(D.expr.trans$CLS)[1] / dim(RRRR$CLS)[1])

plot(ZZZZZZ$edge.wavg)
flaw.ratio4 <- ZZZZZZ$edge.flaws / ZZZZZZ$edge.all.n
plot(flaw.ratio4)

which(flaw.ratio4 <= th4)[1] # 28
which(ZZZZZZ$edge.wavg >= 0.9) # 26
# compressed to k=28

idx <- cutree(H,k=28)
RRRRR <- collapse.cls(CLS, idx, DCOR)
SSSSS <- get.coh.M.score(RRRRR, min.dcor=0.32)

pdf("all.trans.gsplom.k70.k44.k35.k31.k28.d32.apr27.pdf", width=20, height=20)
GGGGG <- splom(RRRRR$CLS, RRRRR$DCOR, grid=F, MIN=0.32, MAX=1, sym=T)
dev.off()

pdf("all.trans.dendro.k70.k44.k35.k31.k28.d32.apr27.pdf", width=20, height=8)
plot(GGGGG$Rhclust)
dev.off()

# ========================================
# compressed to k28. Can we compress even lower?
H <- as.hclust(GGGGG$Rhclust)
CLS <- RRRRR$CLS
DCOR <- RRRRR$DCOR
ZZZZZZZ <- get.compression(CLS,H,DCOR,min.dcor=0.32)
th5 <- 0.05 / (dim(D.expr.trans$CLS)[1] / dim(RRRRR$CLS)[1])

plot(ZZZZZZZ$edge.wavg)
flaw.ratio5 <- ZZZZZZZ$edge.flaws / ZZZZZZZ$edge.all.n
which(flaw.ratio5 <= th5)[1]      # 23
which(ZZZZZZ$edge.wavg >= 0.9)[1] # 26

# compress to k=26
idx <- cutree(H,k=26)
RRRRRR <- collapse.cls(CLS, idx, DCOR)
SSSSSS <- get.coh.M.score(RRRRRR, min.dcor=0.32)

pdf("all.trans.gsplom.k70.k44.k35.k31.k28.k26.d32.apr27.pdf", width=20, height=20)
GGGGGG <- splom(RRRRRR$CLS, RRRRRR$DCOR, grid=F, MIN=0.32, MAX=1, sym=T)
dev.off()

pdf("all.trans.dendro.k70.k44.k35.k31.k28.k26.d32.apr27.pdf", width=20, height=8)
plot(GGGGGG$Rhclust)
dev.off()

# ==================================================
# compressed to k26. compress lower?
H <- as.hclust(GGGGGG$Rhclust)
CLS <- RRRRRR$CLS
DCOR <- RRRRRR$DCOR
ZZZZZZZZ <- get.compression(CLS,H,DCOR,min.dcor=0.32)
th6 <- 0.05 / (dim(D.expr.trans$CLS)[1] / dim(RRRRRR$CLS)[1])

plot(ZZZZZZZZ$edge.wavg)
flaw.ratio6 <- ZZZZZZZZ$edge.flaws / ZZZZZZZZ$edge.all.n
which(flaw.ratio6 <= th6)[1]      # 22
which(ZZZZZZZZ$edge.wavg >= 0.9)[1] # 22

idx <- cutree(H,k=22)
RRRRRRR <- collapse.cls(CLS, idx, DCOR)
SSSSSSS <- get.coh.M.score(RRRRRRR, min.dcor=0.32)

pdf("all.trans.gsplom.k70.k44.k35.k31.k28.k26.k22.d32.apr27.pdf", width=20, height=20)
GGGGGGG <- splom(RRRRRRR$CLS, RRRRRRR$DCOR, grid=F, MIN=0.32, MAX=1, sym=T)
dev.off()

pdf("all.trans.dendro.k70.k44.k35.k31.k28.k26.k22.d32.apr27.pdf", width=20, height=8)
plot(GGGGGGG$Rhclust)
dev.off()


# ==================================================
# compressed to k22. compress lower?
H <- as.hclust(GGGGGGG$Rhclust)
CLS <- RRRRRRR$CLS
DCOR <- RRRRRRR$DCOR
ZY <- get.compression(CLS,H,DCOR,min.dcor=0.32)
th7 <- 0.05 / (dim(D.expr.trans$CLS)[1] / dim(RRRRRRR$CLS)[1])

plot(ZY$edge.wavg)
flaw.ratio7 <- ZY$edge.flaws / ZY$edge.all.n
which(flaw.ratio7 <= th7)[1]      # 19
which(ZY$edge.wavg >= 0.9)[1] # 19

idx <- cutree(H,k=19)
RY <- collapse.cls(CLS, idx, DCOR)
SY <- get.coh.M.score(RY, min.dcor=0.32)

pdf("all.trans.gsplom.k70.k44.k35.k31.k28.k26.k22.k19.d32.apr27.pdf", width=20, height=20)
GY <- splom(RY$CLS, RY$DCOR, grid=F, MIN=0.32, MAX=1, sym=T)
dev.off()

pdf("all.trans.dendro.k70.k44.k35.k31.k28.k26.k22.k19.d32.apr27.pdf", width=20, height=8)
plot(GY$Rhclust)
dev.off()

# ------------------------------------------------------------------------------------------
# compressed to k19. compress lower?
H <- as.hclust(GY$Rhclust)
CLS <- RY$CLS
DCOR <- RY$DCOR
ZYY <- get.compression(CLS,H,DCOR,min.dcor=0.32)
th8 <- 0.05 / (dim(D.expr.trans$CLS)[1] / dim(RY$CLS)[1])

plot(ZYY$edge.wavg)
flaw.ratio8 <- ZYY$edge.flaws / ZYY$edge.all.n
which(flaw.ratio8 <= th8)[1]      # 18
which(ZYY$edge.wavg >= 0.9)[1] # 17

idx <- cutree(H,k=18)
RYY <- collapse.cls(CLS, idx, DCOR)
SYY <- get.coh.M.score(RYY, min.dcor=0.32)

pdf("all.trans.gsplom.k70.k44.k35.k31.k28.k26.k22.k19.k18.d32.apr27.pdf", width=20, height=20)
GYY <- splom(RYY$CLS, RYY$DCOR, grid=F, MIN=0.32, MAX=1, sym=T)
dev.off()

pdf("all.trans.dendro.k70.k44.k35.k31.k28.k26.k22.k19.k18.d32.apr27.pdf", width=20, height=8)
plot(GYY$Rhclust)
dev.off()

# ------------------------------------------------------------------------------------------
# compressed to k18. compress lower?
H <- as.hclust(GYY$Rhclust)
CLS <- RYY$CLS
DCOR <- RYY$DCOR
ZYYY <- get.compression(CLS,H,DCOR,min.dcor=0.32)
th9 <- 0.05 / (dim(D.expr.trans$CLS)[1] / dim(RYY$CLS)[1])

plot(ZYYY$edge.wavg)
flaw.ratio9 <- ZYYY$edge.flaws / ZYYY$edge.all.n
which(flaw.ratio9 <= th9)[1]      # 18
which(ZYY$edge.wavg >= 0.9)[1] # 17
# No. stop.

## ************************************************************
## We have successfully compressed the network to 18 clusters. What is in each of these clusters?
