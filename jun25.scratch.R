# local laptop paths
load("../jun27.GSE2180.gsplom.RData")
source("~/Dropbox/biostat/git_repos/dependency_glyph_splom/lib.R")
# improved greedy GSPLOM dendrogram splitting
# ------------------------------
# 0) filter unrelated elements
# 1) generate symmetric dendrogram
# 2) for k=2 to n
#   find "worst" edge to split
#     break ties
#       a) maximum flaw density (# flaws / # elements)
#       b) most elements
#       c) arbitrarily choose one
#     split both sides, keep split that:
#       a) reduces maximum number of flaws
#       b) results in more specific edges (insig vs sig, directed, sym)
#       c) biggest side
#       d) arbitrary
#   break if no edge meets minimum flaw density or minimum fraction of flawed edges

# for k=2, calculate sum flaws
# str(H) text representation of tree
# attributes(H[[1]]) access named attributes of tree
# dendrapply for applying a function to each node. order.dendrogram and reorder.dendrogram; further, the labels method.

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

make.idx <- function(a,b) {
  as <- paste(as.character(a),sep="",collapse=".")
  bs <- paste(as.character(b),sep="",collapse=".")
  paste(c("",sort(c(as,bs)),""),sep="",collapse="|")
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

# initialize function parameters.
DC <- DCOR.TF
BC <- BOOL.TF
H <- R.GSE2180.TF$Rhclust
sig <- 0.36
# initialize sub-divisions as left and right branches of dendrogram root.
TREE <- list(c(1), c(2))
# initialize flaw counts as the single edge between the left and right sub-branches
FLAWS <- list(); SIZES <- list(); DENSITIES <- list()
i <- make.idx(c(1),c(2))
FLAWS[[i]] <- get.num.flaws(DC,BC,H,c(1),c(2),sig)
SIZES[[i]] <- count.members(H,c(1),c(2))
DENSITIES[[i]] <- FLAWS[[i]]/SIZES[[i]]

# 1 find max dense edge
max.dense.i <- which.max(DENSITIES)
# 2 split rows and columns; find which maximially reduces number of flaws



branch1<-c(1)
branch2<-c(2)
qq.1 <- order.dendrogram(H[[branch1]])
qq.2 <- order.dendrogram(H[[branch2]])
BS <- BC[qq.1,qq.2]
DS <- DC[qq.1,qq.2]
sig.mask <- DS>=sig

summary(as.factor(BC[qq.2,qq.1]))
#    0     1     2     3     4     5     6     7 
#19809  6915   493   807  9083  2762    18   127 
summary(as.factor(BC[qq.1,qq.2]))
#    0     1     2     3     4     5     6     7 
#19809   807   493  6915  9083  2762    18   127
summary(as.factor(BS))
##     0     1     2     3     4     5     6     7 
## 19809  6915   493   807  9083  2762    18   127 
summary(as.factor(BS[sig.mask]))
##    0    1    2    3    4    5    6    7 
## 9169 4796  493  692 2794 1360   18   65 
count.mix.sign(BS[sig.mask])
#[1] 1443
count.mix.pos.dir(BS[sig.mask])
#[1] 692
count.mix.neg.dir(BS[sig.mask])
#[1] 65

FLAWS[[i]] <- get.num.flaws(DC,BC,H,c(1),c(2),sig)
# [1] 1881
SIZES[[i]] <- count.members(H,c(1),c(2))
# 40014
DENSITIES[[i]] <- FLAWS[[i]]/SIZES[[i]]
# $`|1|2|`
# [1] 0.04700855
