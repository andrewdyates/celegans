# Filter probes to genes. Re-annotate.
# ------------------------------
library(Biobase)
load("../celegans.apr8.expr.RData")
load("../apr9.dcor.cls.RData")
# no _x_ or missing gene probes
probe.to.sym <- read.table("../GSE2180_GPL200.probe2genesym.tab", sep="\t", stringsAsFactors=F, row.names=1)
trans.factors <- unique(as.vector(unlist(read.table("celegans.transcription.factors.txt", stringsAsFactors=F))))
gold.genes <- c('pal-1', 'tbx-8', 'tbx-9', 'elt-1', 'hnd-1', 'scrt-1', 'cwn-1', 'unc-120', 'hlh-1', 'nob-1', 'elt-3', 'nhr-25', 'mab-21', 'lin-26', 'vab-7')

# Super-sets of the probe sets that we want
# E.expr, E.expr.trans, E.gold

# Assign gene symbol to featureData
qq <- match(rownames(probe.to.sym), rownames(exprs(E.expr)))
Eg.expr <- E.expr[qq[!is.na(qq)], ]
qq2 <- match(rownames(exprs(Eg.expr)), rownames(probe.to.sym))
syms <- probe.to.sym[qq2[!is.na(qq2)],]
featureData(Eg.expr)$genesyms <- syms

# choose highest expressed probe per gene
# choose highest mean, non _s_ probe. If all _s_, choose highest _s_.
means <- rowMeans(exprs(Eg.expr))
mean.exps <- split(means,featureData(Eg.expr)$genesyms)

choose.best.probe <- function(s) {
  probes <- sort(s, decreasing=T)
  r <- names(probes)[1]
  # if r is an _s_ probe, choose the first non _s_ probe if it exists
  qq <- grep("_s_", names(probes), invert=T)
  if(length(grep("_s_", r)) && length(qq)) 
    r <- names(probes)[qq[1]]
  r
}
best.gene.probes <- lapply(mean.exps, function(s) choose.best.probe(s))
# [1] 5336 representative, expressed probes
qq <- rownames(exprs(Eg.expr)) %in% best.gene.probes
Eg.expr <- Eg.expr[qq,]
v <- sapply(featureData(Eg.expr)$genesyms, function(s) strsplit(s,";"))
featureData(Eg.expr)$genesyms.v <- v

# ----------------------------------------
# Choose Best Probes per Gene list
# ----------------------------------------
get.probes <- function(gene.list) {
  probe.i <- sapply(gene.list, function(g) {
    which(sapply(featureData(Eg.expr)$genesyms.v, function(v) any(g %in% v)))
  })
  lapply(probe.i, function(z) {
    if(length(z) == 1)
      as.vector(unlist(z))
    if(length(grep(";", names(z))) == length(z))
      # no good choice, just choose the top one
      as.vector(unlist(z))[1]
    else
      # choose name without multiple genes
      z[[grep(";", names(z), invert=T)[1]]]
  })
}

gold.i <- get.probes(gold.genes)
trans.i <- get.probes(trans.factors)
length(trans.factors)
#[1] 609
sum(!is.na(trans.i))
#[1] 191

qq <- gold.i[!is.na(gold.i)]
Eg.expr.gold <- Eg.expr[unlist(qq),]
featureData(Eg.expr.gold)$sym <- names(gold.i)[!is.na(gold.i)]

qq <- trans.i[!is.na(trans.i)]
Eg.expr.trans <- Eg.expr[unlist(qq),]
featureData(Eg.expr.trans)$sym <- names(trans.i)[!is.na(trans.i)]

save(Eg.expr.gold, Eg.expr.trans, Eg.expr, file="../apr16.genelevel.exprs.RData")
