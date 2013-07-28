# local laptop paths
load("../jun27.GSE2180.gsplom.RData")
source("~/Dropbox/biostat/git_repos/dependency_glyph_splom/lib.R")
# improved greedy GSPLOM dendrogram splitting
# ------------------------------
# 0) filter unrelated elements
# 1) generate symmetric dendrogram
# 2) from k=2, split along the dendrogram branch that resolves the most contradictions, else split biggest

# Verify that all labels are aligned as expected
stopifnot(sum(duplicated(rownames(DCOR.TF)))==0)
stopifnot(all(rownames(DCOR.TF) %in% labels(R.GSE2180.TF$Rhclust)))
stopifnot(all(labels(R.GSE2180.TF$Rhclust) %in% rownames(DCOR.TF)))
stopifnot(all(labels(R.GSE2180.TF$Rhclust) == rownames(DCOR.TF)[order.dendrogram(R.GSE2180.TF$Rhclust)]))
stopifnot(all(!order.dendrogram(R.GSE2180.TF$Rhclust[[c(1)]]) %in% order.dendrogram(R.GSE2180.TF$Rhclust[[c(2)]])))
stopifnot(all(labels(R.GSE2180.TF$Rhclust[[1]]) == rownames(DCOR.TF)[order.dendrogram(R.GSE2180.TF$Rhclust[[1]])]))
stopifnot(all(labels(R.GSE2180.TF$Rhclust[[2]]) == rownames(DCOR.TF)[order.dendrogram(R.GSE2180.TF$Rhclust[[2]])]))

count.flaws <- function(CLS, DCOR, sig=0) {
  sig.mask <- DCOR>=sig
  flaws <- count.mix.sign(CLS[sig.mask]) + count.mix.pos.dir(CLS[sig.mask]) + count.mix.neg.dir(CLS[sig.mask])
  flaws
}

count.members <- function(H,b1,b2) {
  attributes(H[[b1]])$members * attributes(H[[b2]])$members
}

make.idx.single <- function(a) {
  as <- paste(as.character(a),sep="",collapse=".")
}
make.idx <- function(a,b) {
  as <- paste(as.character(a),sep="",collapse=".")
  bs <- paste(as.character(b),sep="",collapse=".")
  paste(c("",sort(c(as,bs)),""),sep="",collapse="|")
}
rev.idx.single <- function(i) {
  as.numeric(unlist(strsplit(i,'.',fixed=TRUE)))
}
rev.idx <- function(i) {
  z <- strsplit(substr(i,2,nchar(i)),'|',fixed=TRUE)[[1]]
  i1 <- as.numeric(unlist(strsplit(z[1],'.',fixed=TRUE)))
  i2 <- as.numeric(unlist(strsplit(z[2],'.',fixed=TRUE)))
  list(i1,i2)
}

get.num.flaws <- function(DC,BC,H,branch1,branch2,sig=0) {
  # select row/col numbers from branch1, branch2 in dendrogram H
  qq.1 <- order.dendrogram(H[[branch1]])
  qq.2 <- order.dendrogram(H[[branch2]])
  # from subselection of these row/col numbers, return number of flaws
  count.flaws(BC[qq.1,qq.2], DC[qq.1,qq.2], sig=sig)
}

# For each sub-tree, track how many flaws would be resolved
# if I split on sub-tree A given sub-trees A and B...
test.splits.hasflaw <- function(TREES,DC,BC,H,sig) {
  R <- list()
  for (As in TREES) {
    A <- rev.idx.single(As)
    has.flaws <- as.numeric(get.num.flaws(DC,BC,H,c(A,1),c(A,2),sig)>0)
    for (Bs in TREES) {
      if (As == Bs) next
      B <- rev.idx.single(Bs)
      has.flaws <- has.flaws + as.numeric(get.num.flaws(DC,BC,H,c(A,1),B,sig)>0)
      has.flaws <- has.flaws + as.numeric(get.num.flaws(DC,BC,H,c(A,2),B,sig)>0)
    }
    R[[make.idx.single(A)]] <- has.flaws
  }
  R
}
test.splits.flawdense <- function(TREES,DC,BC,H,sig) {
  R <- list()
  for (As in TREES) {
    A <- rev.idx.single(As)
    sum.wfd <- 0
    for (Bs in TREES) {
      if (As == Bs) next
      B <- rev.idx.single(Bs)
      f <- get.num.flaws(DC,BC,H,A,B,sig)
      n <- count.members(H,A,B)
      f1 <- get.num.flaws(DC,BC,H,c(A,1),B,sig)
      n1 <- count.members(H,c(A,1),B)
      f2 <- get.num.flaws(DC,BC,H,c(A,2),B,sig)
      n2 <- count.members(H,c(A,2),B)
      fd <- f/n - (n1/n)*(f1/n1)+(n2/n)*(f2/n2)
      wfd <- (attributes(H[[B]])$members/n) * fd
      sum.wfd <- sum.wfd + wfd
    }
    R[[make.idx.single(A)]] <- sum.wfd
  }
  R
}

test.splits.flawdense2 <- function(TREES,DC,BC,H,sig) {
  R <- list()
  for (As in TREES) {
    A <- rev.idx.single(As)
    sum.wfd <- 0
    for (Bs in TREES) {
      if (As == Bs) next
      B <- rev.idx.single(Bs)
      f <- get.num.flaws(DC,BC,H,A,B,sig)
      n <- count.members(H,A,B)
      wfd <- (attributes(H[[B]])$members/n) * (f/n)
      sum.wfd <- sum.wfd + wfd
    }
    R[[make.idx.single(A)]] <- sum.wfd
  }
  R
}

branch.sizes <- function(TREES, H) {
  R <- list()
  for (A in TREES) {
    A <- rev.idx.single(A)
    R[[make.idx.single(A)]] <- attributes(H[[A]])$members
  }
  R
}

# initialize function parameters.
DC <- DCOR.TF
BC <- BOOL.TF
H <- R.GSE2180.TF$Rhclust
sig <- 0.36

# always start with a split root
TREES <- c("1","2")

#R.SPLITS <- test.splits.hasflaw(TREES, DC, BC, H, sig)
R.SIZES <- branch.sizes(TREES,H)

iv <- which.min(R.SPLITS)
# if multiple, choose largest
if (sum(R.SPLITS == R.SPLITS[[iv]]) > 1) {
  iv <- which.max(R.SIZES[R.SPLITS == R.SPLITS[[iv]]])
}

# split selection
niv <- names(iv)
TREES <- TREES[TREES!=niv]
TREES <- c(TREES, paste0(niv,".1"), paste0(niv,".2"))


# ------------------------------
TREES <- c("1","2")
R.FLAWDENSE <- test.splits.flawdense(TREES, DC, BC, H, sig)
R.SIZES <- branch.sizes(TREES,H)
iv <- which.max(R.FLAWDENSE)
# if multiple, choose largest
if (sum(R.FLAWDENSE == R.FLAWDENSE[[iv]]) > 1) {
  iv <- which.max(R.SIZES[R.FLAWDENSE == R.FLAWDENSE[[iv]]])
}

# split selection
niv <- names(iv)
TREES <- TREES[TREES!=niv]
TREES <- c(TREES, paste0(niv,".1"), paste0(niv,".2"))
R.FLAWDENSE
TREES

# This algorithm doesn't work very well...
