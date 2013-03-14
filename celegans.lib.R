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
all.pairs.dcov <- function(M) {
  n <- nrow(M)
  D <- outer(1:n,1:n, FUN = Vectorize(function(i,j) dcov(M[i,],M[j,])))
  rownames(D) <- rownames(M)
  colnames(D) <- rownames(M)
  D
}

all.steps <- function(M, do.plot=TRUE) {
  STEPS <- apply(M, 1, fit.upstep)
  if(do.plot) {
    for (i in 1:dim(M)[1]) {
      title=rownames(M)[i]
      plot.stepfit(STEPS[[i]], v=M[i,], add.mean.median=T, main=paste(title, "Stepfit"))
      plot.sse(STEPS[[i]], add.mean.median=T, main=paste(title, "SSE"))
    }
  }
  STEPS
}

all.pairs.cls <- function(M, steps, b=0.5, do.plot=TRUE) {
  n <- dim(M)[1]
  CLS <- mat.or.vec(n,n)
  for(i in 1:n) { # row
    for(j in 1:n) { # col
      y <- M[i,]
      y.th <- steps[[i]]$th
      x <- M[j,]
      x.th <- steps[[j]]$th
      x.title=rownames(M)[j]
      y.title=rownames(M)[i]
      RR <- cls.pair(x,y,x.th,y.th, b.x=b, b.y=b, do.plot=(i<j)&&do.plot, xlab=x.title, ylab=y.title)
      CLS[i,j] <- cls.to.enum(RR$CLS)
    }
  }
  rownames(CLS) <- rownames(M)
  colnames(CLS) <- rownames(M)
  CLS
}


MIN <- 0
MAX <- 1
cols <- brewer.pal(8,"RdYlBu")
heatmap_breaks <- seq(MIN,MAX,0.01)
heatmap_cols <- rev(colorRampPalette(cols)(length(heatmap_breaks)-1))

heatmap_breaks.near0 <- seq(-0.4,0.4,0.01)
heatmap_cols.near0 <- rev(colorRampPalette(cols)(length(heatmap_breaks.near0)-1))

gold.genes <- c("pal-1", "tbx-8", "tbx-9", "elt-1", "cwn-1", "hnd-1", "scrt-1", "elt-3", "nhr-25", "nob-1", "vab-7", "unc-120", "hlh-1", "mab-21", "lin-26")
known.noexpress <- c("193640_s_at", "175709_at", "193846_at")


