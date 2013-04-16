# Filter probes to genes. Re-annotate.
# ------------------------------
library(Biobase)
load("../celegans.apr8.expr.RData")
load("../apr9.dcor.cls.RData")
probe.to.sym <- read.table("../GSE2180_GPL200.probe2genesym.tab", sep="\t", stringsAsFactors=F, row.names=1)
trans.factors <- unique(as.vector(unlist(read.table("celegans.transcription.factors.txt", stringsAsFactors=F))))

# Super-sets of the probe sets that we want
# E.expr, E.expr.trans, E.gold
# Assign gene symbol to featureData
qq <- match(rownames(exprs(E.expr)), rownames(probe.to.sym))



