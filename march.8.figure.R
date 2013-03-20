library(sva)
library(lumi) # for nice figures, plus Biobase
library(Biobase)
library(energy)

gold.genes <- c("pal-1", "tbx-8", "tbx-9", "elt-1", "cwn-1", "hnd-1", "scrt-1", "elt-3", "nhr-25", "nob-1", "vab-7", "unc-120", "hlh-1", "mab-21", "lin-26")

GSE9665.GEO <- read.table("../GSE9665_GPL200.seriesmatrix.tab", sep="\t", header=TRUE, row.names=1, comment="")
GSE2180.GEO <- read.table("../GSE2180_GPL200.seriesmatrix.tab", sep="\t", header=TRUE, row.names=1, comment="")
GPL200 <- read.table("../GPL200-2880.txt", sep="\t", quote="", comment="", header=TRUE, row.names=1)

## Load GSE2180 phenotype data.
## ------------------------------
AttrT <- read.table("../GSE2180_GPL200.samples.tab", sep="\t", header=TRUE, row.names=1)
AttrT <- as.data.frame(t(AttrT))
AttrT$time <- as.numeric(as.character(AttrT$"n:Sample_time"))
AttrT$genotype <- AttrT$"n:Sample_genotype"
GSE2180.Samples <- AttrT

## Load GSE9665 phenotype data.
## ------------------------------
AttrT <- read.table("../GSE9665_GPL200.samples.tab", sep="\t", header=TRUE, row.names=1)
AttrT <- as.data.frame(t(AttrT))
AttrT$knockout <- sub("mex-3 (\\w+-[0-9,]+)? ?(\\d+) (\\w).*", AttrT$title, replacement="\\1")
AttrT$knockout[AttrT$knockout==""] <- "WT"
AttrT$knockout <- as.factor(AttrT$knockout)
AttrT$num <- as.factor(sub("mex-3 (\\w+-[0-9,]+)? ?(\\d+) (\\w).*", AttrT$title, replacement="\\2"))
AttrT$primer <- as.factor(sub("mex-3 (\\w+-[0-9,]+)? ?(\\d+) (\\w).*", AttrT$title, replacement="\\3"))
GSE9665.Samples <- AttrT

#
qq.c <- match(colnames(GSE2180.GEO), rownames(GSE2180.Samples))
qq.r <- match(rownames(GSE2180.GEO), rownames(GPL200))
GSE2180 <- ExpressionSet(
  as.matrix(GSE2180.GEO),
  AnnotatedDataFrame(GSE2180.Samples[qq.c,]),
  AnnotatedDataFrame(GPL200[qq.r,]))

qq.c <- match(colnames(GSE9665.GEO), rownames(GSE9665.Samples))
qq.r <- match(rownames(GSE9665.GEO), rownames(GPL200))
GSE9665 <- ExpressionSet(
  as.matrix(GSE9665.GEO),
  AnnotatedDataFrame(GSE9665.Samples[qq.c,]),
  AnnotatedDataFrame(GPL200[qq.r,]))

save(GSE2180, GSE9665, file="../geo.GSE2180.GSE9665.RData")
## ========================================


## ----------------------------------------
## GOLD STANDARD NETWORKS
## ----------------------------------------
# GSE2180 Time Series
row.nums <- sapply(gold.genes, function(s) grep(paste0("(",s," |",s,"$)"), featureData(GSE2180)$Gene.Symbol))
## visually inspect, get correct probe IDs
sapply(row.nums, function(q) featureData(GSE2180)$Gene.Symbol[q]) ## OK
GSE2180.GOLD <- GSE2180[unlist(row.nums),]

pdf("gse2180.geo.gold.pairs.pdf", width=30, height=30)
pairs(t(exprs(GSE2180.GOLD)), labels=names(unlist(row.nums)), main="GSE2180 C.Elegans GSN Gene Names (with multi-probe enumerations)")
pairs(t(exprs(GSE2180.GOLD)), main="GSE2180 C.Elegans GSN Probes")
dev.off()

# ----------
# GSE9665 Knockouts
row.nums.2 <- sapply(gold.genes, function(s) grep(paste0("(",s," |",s,"$)"), featureData(GSE9665)$Gene.Symbol))
## visually inspect, get correct probe IDs
sapply(row.nums.2, function(q) featureData(GSE9665)$Gene.Symbol[q]) ## OK
GSE9665.GOLD <- GSE9665[unlist(row.nums.2),]

pdf("gse9665.geo.gold.pairs.pdf", width=30, height=30)
pairs(t(exprs(GSE9665.GOLD)), labels=names(unlist(row.nums)), main="GSE9665 C.Elegans GSN Gene Names (with multi-probe enumerations)")
pairs(t(exprs(GSE9665.GOLD)), main="GSE9665 C.Elegans GSN Probes")
dev.off()
## ========================================

## VISUALIZE SPLOMS
source("~/Dropbox/biostat/git_repos/dependency_glyph_splom/lib.R")
source("~/Dropbox/biostat/git_repos/boolean_implication_fit/bool.R")
source("~/Dropbox/biostat/git_repos/boolean_implication_fit/step.up.R")
source("~/Dropbox/biostat/git_repos/boolean_implication_fit/run.all.bool.R")

all.pairs.dcor <- function(M) {
  n <- nrow(M)
  outer(1:n,1:n, FUN = Vectorize(function(i,j) dcor(M[i,],M[j,])))
}

