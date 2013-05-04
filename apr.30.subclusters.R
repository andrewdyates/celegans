library(Biobase)
load("../trans.genes.apr17.gsplom.RData")
load("../apr16.genelevel.depM.RData")
load("../apr16.genelevel.exprs.RData")
source("~/Dropbox/biostat/git_repos/dependency_glyph_splom/lib.R")
source("collapse.weak.R")
source("all.pairs.weak.R")
source("compress.classes.R")
load(file="../altrans.weak.RData")

## Plot compression measures
load("quant.all.trans.k.RData")

flaw.fract <- ZZ$edge.flaws / ZZ$edge.all.n
pdf("all.trans.flaw.fract.pdf", width=5, height=5)
plot(flaw.fract, main="Edge Flaw Fraction Per Number of Clusters", ylab="Flaw Fraction", xlab="Number of Clusters")
abline(h=0.05, col="#0093c3", lty=2, lwd=2)
abline(v=61, col="#0093c3", lty=1, lwd=2)
dev.off()

pdf("number.dcor0.32.edges.pdf", width=5, height=5)
plot(ZZ$edge.n, ylab="Number Edges >= dCor 0.32", xlab="Number of Clusters")
abline(v=61, col="#0093c3", lty=1, lwd=2)
abline(h=ZZ$edge.n[61], col="#0093c3", lty=2, lwd=2)
dev.off()

max(ZZ$edge.n) # 12225
ZZ$edge.n[61] # 1108
pdf("fraction.significant.edges.dcor0.32.pdf", width=5, height=5)
plot(ZZ$edge.n/ZZ$edge.all.n, main="Fraction of significant edges at dCor 0.32", xlab="Number of Clusters", ylab="Fraction") # dCor 0.32
abline(v=61, col="#0093c3", lty=1, lwd=2)
abline(h=(ZZ$edge.n/ZZ$edge.all.n)[61], col="#0093c3", lty=2, lwd=2)
dev.off()

H <- as.hclust(TRANS.R$Rhclust)
CLS <- D.expr.trans$CLS
DCOR <- D.expr.trans$DCOR
idx <- cutree(H,k=61)
R61 <- collapse.cls(CLS, idx, DCOR)
S61 <- get.coh.M.score(R61, min.dcor=0.32)

# 503 edges to clusters
summary(c(R61$DCOR.EXPAND[upper.tri(DCOR)]))

qq <- is.na(R61$DCOR.EXPAND[upper.tri(DCOR)])
sum(DCOR[upper.tri(DCOR)][qq] >= 0.32) # 503, all are
sum(CLS[upper.tri(CLS)][qq] %in% c(1,2,3,4))
# how many edges total?
zz <- ((DCOR >= 0.32) & (CLS %in% c(1,2,3,4)))[upper.tri(DCOR)]

# how many edges in r61 compression?
d <- R61$DCOR[upper.tri(R61$DCOR)] >= 0.32
c <- as.factor(R61$CLS[upper.tri(R61$CLS)])
d & c ==1

# get class / edge change
# how many edges in full are no longer edges in k61 compression?
d <- DCOR[upper.tri(DCOR)]
d61 <- R61$DCOR.EXPAND[upper.tri(R61$DCOR.EXPAND)]
d61.na <- is.na(d61)
df <- d[!d61.na]
d61f <- d61[!d61.na]
clsf <- CLS[upper.tri(CLS)][!d61.na]
cls61f <- R61$CLS.EXPAND[upper.tri(CLS)][!d61.na]

sum(df >= 0.32 & d61f < 0.32) # 1029
sum(df < 0.32 & d61f >= 0.32) # 1774
sum(df >= 0.32 & d61f >= 0.32) # 10693
sum(df < 0.32 & d61f < 0.32)  # 4146
# total: 17642
# 17642 + 503 == 1814

lost.cls.t <- as.factor(clsf[df >= 0.32 & d61f < 0.32])
lost.cls.61 <- as.factor(cls61f[df >= 0.32 & d61f < 0.32])

gain.cls.t <- as.factor(clsf[df < 0.32 & d61f >= 0.32])
gain.cls.61 <- as.factor(cls61f[df < 0.32 & d61f >= 0.32])

# of the edges that are significant in in both, did their classes change?
same.cls.t <- as.factor(clsf[df >= 0.32 & d61f >= 0.32])
same.cls.61 <- as.factor(cls61f[df >= 0.32 & d61f >= 0.32])

# row: original class. col: new class
CNT <- matrix(0, nrow=8, ncol=8)
c.alls <- c(0,1,2,3,4,5,6,7)
c.61s <- c(8,1,2,3,4,5,6,7)
for (i.all in 1:8) {
  for (i.61 in 1:8) {
    c.all <- c.alls[i.all]
    c.61 <- c.61s[i.61]
    CNT[i.all,i.61] = sum((same.cls.t==c.all) & (same.cls.61==c.61))
  }
}



# difference between compressed and actual dCOR value
D <- (R61$DCOR.EXPAND - DCOR)[upper.tri(DCOR)]
hist(D)
sd(D,na.rm=T)
#[1] 0.08248709
D[is.na(D)] <- 0
> sd(D)
#[1] 0.08133567
pdf("../dcor_diff_hist.pdf", width=6, height=6)
hist(D)
dev.off()

summary(D)
length(D) #18145 edges
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#-0.3880 -0.0489  0.0014  0.0000  0.0500  0.3546     503
## what about significant edges only?
sum(DCOR >= 0.32)
#[1] 24641
sum(DCOR >= 0)
#[1] 36481
D.sig <- (R61$DCOR.EXPAND - DCOR)[upper.tri(DCOR) & DCOR >= 0.32]
length(D.sig)
# 12225
summary(D.sig)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#-0.3880 -0.0737 -0.0216 -0.0216  0.0291  0.3546     503 
#> hist(D.sig)
sd(D.sig,na.rm=T)
#[1] 0.08495793
## number of edges from significant to non-significant
D.sig2non <- DCOR >= 0.32 & R61$DCOR.EXPAND < 0.32
sum(D.sig2non[upper.tri(D.sig2non)], na.rm=T)
#[1] 1029
sum(upper.tri(D.sig2non)) # total number of edges
# [1] 18145
1029/18145
#[1] 0.05670984 fraction of significant edges now below threshold
D.non2sig <- DCOR < 0.32 & R61$DCOR.EXPAND >= 0.32
sum(D.non2sig[upper.tri(D.non2sig)], na.rm=T)
#[1] 1774
1774/18145
#[1] 0.09776798
D.nonsig <- (R61$DCOR.EXPAND - DCOR)[upper.tri(DCOR) & DCOR < 0.32]
length(D.nonsig)
# 5920 (note: 12225 sig edges, 48.4% edges significant); 5920 + 12225 == 18145
# @ 0.5 threshold, there are only 36 implied relationships below original threshold 
D.non2sig <- DCOR < 0.32 & R61$DCOR.EXPAND >= 0.5
sum(D.non2sig[upper.tri(D.non2sig)], na.rm=T)
#> [1] 36
#> 36 / 18145
#[1] 0.001984018
# total possible implied inter-cluster edges
sum(!is.na(R61$DCOR.EXPAND[upper.tri(DCOR)]))
#[1] 17642
17642/18145
#[1] 0.9722789 @ k61


## for significant edges, how frequently does an edge change?
summary(as.factor(c(R61$CLS.EXPAND)))
 ##    1     2     3     4     5     7     8  NA's 
 ## 4941   366  4941 15788  8408    16   824  1197 
summary(as.factor(c(CLS))) # note: NA is 0 in CLS, 8 in EXPAND
##    0     1     2     3     4     5     6     7 
## 1326  5355  1503  5355 14308  8548    28    58 

# ====
# edges that change at original 0.32?
CLS.q <- upper.tri(CLS) & (DCOR >= 0.32)

# how to display this?
q1 <- CLS[CLS.q] == 1
sum1 <- summary(as.factor(R61$CLS.EXPAND[CLS.q][q1]))
q2 <- CLS[CLS.q] == 2
sum2 <- summary(as.factor(R61$CLS.EXPAND[CLS.q][q2]))
q3 <- CLS[CLS.q] == 3
sum3 <- summary(as.factor(R61$CLS.EXPAND[CLS.q][q3]))
q4 <- CLS[CLS.q] == 4
sum4 <- summary(as.factor(R61$CLS.EXPAND[CLS.q][q4]))
qnon <- CLS[CLS.q] %in% c(5,6,7,0)
sumnon <- summary(as.factor(R61$CLS.EXPAND[CLS.q][qnon]))

# ====
# edges that change at summarized 0.5?
CLS.q <- upper.tri(CLS) & (R61$DCOR.EXPAND >= 0.5) & !is.na(R61$DCOR.EXPAND)

# how to display this?
q1 <- CLS[CLS.q] == 1
sum1 <- summary(as.factor(R61$CLS.EXPAND[CLS.q][q1]))
q2 <- CLS[CLS.q] == 2
sum2 <- summary(as.factor(R61$CLS.EXPAND[CLS.q][q2]))
q3 <- CLS[CLS.q] == 3
sum3 <- summary(as.factor(R61$CLS.EXPAND[CLS.q][q3]))
q4 <- CLS[CLS.q] == 4
sum4 <- summary(as.factor(R61$CLS.EXPAND[CLS.q][q4]))
qnon <- CLS[CLS.q] %in% c(5,6,7,0)
sumnon <- summary(as.factor(R61$CLS.EXPAND[CLS.q][qnon]))


#flaw.fract[61]
#[1] 0.04754098
## Get edge transition matrix
