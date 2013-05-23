# 0: no class; 1: and; 2: rn4c (row necessary for col); 3: cn4r (col necessary for row); 4: xor; 5: mix, 6: no class

collapse.weak <- function(WEAK, idx) {
  idx <- as.factor(idx)
  n <- length(levels(idx))
  COMP.W <- matrix(0, nrow=n, ncol=n)
  rownames(COMP.W) <- 1:n
  colnames(COMP.W) <- 1:n
  for (i in 1:length(levels(idx))) {
    for (j in 1:length(levels(idx))) {
      iv <- idx == levels(idx)[i]
      jv <- idx == levels(idx)[j]
      W <- WEAK[iv,jv]
      w <- choose.weak(W, sym=i==j)
      COMP.W[i,j] <- w
    }
  }
  COMP.W
}

choose.weak <- function(WEAK, sym=F) {
  cnt <- weak.counts(WEAK, sym=sym)
  n <- sum(cnt)
  max.cnt <- max(cnt)
  max.cnt.i <- cnt == max.cnt
  if (max.cnt < n/2) {
    return(5)
  } else if (max.cnt.i[2]) {
    if (cnt[3] == 0 && cnt[4] == 0)
      return(2)
    else
      return(5)
  } else if (max.cnt.i[3]) {
    if (cnt[2] == 0 && cnt[4] == 0)
      return(3)
    else
      return(5)
  } else if (max.cnt.i[4]) {
    if(cnt[2] == 0 && cnt[3] == 0 && cnt[1] == 0)
      return(4)
    else
      return(5)
  } else if (max.cnt.i[1]) {
    if(cnt[4] == 0)
      return(1)
    else
      return(5)
  } else {
    return(5)
  }
  print("WARNING: edge case in choose.weak.")
  5
}
      
weak.counts <- function(WEAK, sym=F) {
  WEAK[WEAK==0] <- 6
  if(sym) {
    stopifnot(dim(WEAK)[1] == dim(WEAK)[2])
    W <- WEAK[upper.tri(WEAK)]
  } else {
    W <- c(WEAK)
  }
  weak.counts <- sapply(c(1:6), function(i) sum(W==i))
  # account for arbitrary symmetry across diagonal
  if (sym) {
    weak.counts[2] = weak.counts[2] + weak.counts[3]
    weak.counts[3] = 0
  }
  weak.counts
}
