load("../jul1.NAfiltered.tf.RData")
source("~/Dropbox/biostat/git_repos/dependency_glyph_splom/lib.R")
pdf("~/Desktop/GSE2180.gsplom.tf.filt.crits.jul15.k420.pdf", width=100, height=100)
R.GSE2180.F.TF <- splom(BOOL.TF.F, DCOR.TF.F, useRaster=T, sym=T, high.sat=T, MIN=0.36, MAX=1, row.cls.dist=BOOL.TF.F.D)
dev.off()

qq <- order.dendrogram(R.GSE2180.F.TF$Rhclust)
gsplom.names <- rownames(DCOR.TF.F)[qq]

writeLines(gsplom.names, "../sep30.GSPLOM.names.txt")
save(R.GSE2180.F.TF, gsplom.names, file="../sep30.R.GSE2180.F.TF")
