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
    u <- DCOR # mean is value itself
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
  R$CLS.EXPAND <- matrix(NA, nrow=dim(CLS)[1], ncol=dim(CLS)[2])
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
    R$DCOR.EXPAND <- matrix(NA, nrow=dim(DCOR)[1], ncol=dim(DCOR)[2])
  } else {
    R$DCOR <- NULL
  }
  # compute overall coherence
  for (i in 1:length(levels(idx))) {
    for (j in 1:length(levels(idx))) {
      iv <- idx == levels(idx)[i]
      jv <- idx == levels(idx)[j]
      C <- CLS[iv,jv]
      if (!is.null(DCOR)) {
        R$DCOR[i,j] <- get.mean.dcor(DCOR[iv,jv], sym=i==j)
        if (i!=j) {
          R$DCOR.EXPAND[iv,jv] <- R$DCOR[i,j]
        }
      }
      r <- choose.coh.cls(get.coh.vec(C, sym=i==j))
      if (is.null(r$cls)) {
        print(C)
        print(r)
        stop("Null class.")
      }
      R$CLS[i,j] <- r$cls
      if (i!=j) {
        R$CLS.EXPAND[iv,jv] <- r$cls
      }
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
  R$pos.edge.n <- sum(COLLAPSED$CLS[tri] %in% c(1,2,3,4))
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
    R$pos.edge.n[k] <- S$pos.edge.n
    R$edge.all.n[k] <- S$edge.all.n
    R$clust.n[k] <- S$clust.n
    R$clust.all.n[k] <- S$clust.all.n
    # how many
    if (!is.null(S$dcor.clust))
      R$dcor.clust[k] <- S$dcor.clust
  }
  R
}
