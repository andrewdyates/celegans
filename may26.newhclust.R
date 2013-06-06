library(Biobase)
load("../trans.genes.apr17.gsplom.RData")
load("../apr16.genelevel.depM.RData")
load("../apr16.genelevel.exprs.RData")
source("~/Dropbox/biostat/git_repos/dependency_glyph_splom/lib.R")


system.time(D.cls.r <- gen.glyph.dist.m(D.expr.trans$CLS))
# 85.153   0.186  86.692 
system.time(D.dcor.r <- dist(D.expr.trans$DCOR))
#  0.013   0.000   0.014 
D.row <- D.dcor.r*2+sqrt(D.cls.r)
system.time(H <- as.dendrogram(hclust(D.row, method="average")))
#  0.014   0.000   0.014 
library("flashClust")
system.time(H <- as.dendrogram(hclust(D.row, method="average")))
#  0.015   0.000   0.015 

