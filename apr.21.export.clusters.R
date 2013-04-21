library(Biobase)
load("../apr16.genelevel.depM.RData")
load("../apr16.genelevel.exprs.RData")
source("~/Dropbox/biostat/git_repos/dependency_glyph_splom/lib.R")


pdf("/dev/null")
GOLD.R <- splom(D.expr.gold$CLS, D.expr.gold$DCOR, MIN=0.3, MAX=1)
dev.off()

H <- as.hclust(GOLD.R$Rhclust)
H.clusts <- cutree(H, k=7)

write.table(H.clusts, file="gold.k7.hclust.csv", sep=",", col.name=F)
