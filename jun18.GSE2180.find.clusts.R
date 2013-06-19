source("~/code/dependency_glyph_splom/lib.R")
load("../jun13.GSE2180.gsplom.RData")
GSE2180.Z <- get.compression(BOOL.TF, as.hclust(R.GSE2180.TF$Rhclust), DCOR.TF, min.dcor=0.357, max.k=300)
save(GSE2180.Z, file="../jun18.GSE2180.Z.compressions.pdf")

