library("Biobase")
library("RColorBrewer")
library("gplots")
library("dynamicTreeCut")
load("../celegans.apr8.expr.RData")
load("../apr9.dcor.cls.RData")
source("~/pymod/dependency_glyph_splom/lib.R")

# color gold standard probes
phaseI <- c("pal-1 174043_at", "pal-1 193341_at", "pal-1 193342_s_at")
phaseII.out <- c("tbx-9 190539_s_at", "scrt-1 192307_at")
phaseII <- c("cwn-1 188486_s_at", "elt-1 192655_s_at", "hnd-1 192707_at", "tbx-8 190477_at")
phaseIII <- c("elt-3 175801_at", "nob-1 175771_at", "nob-1 191468_s_at", "nhr-25 175663_at", "unc-120 188706_at", "hlh-1 193759_at")
phaseIV <- c("mab-21 188239_s_at")
phaseNA <- c("lin-26;lir-1 174037_at", "lin-26;lir-1 190309_at")
cols <- brewer.pal(6,"Set1")
classCols <- rep("#ffffff", dim(TRANS.CLS)[1])
pns <- rownames(TRANS.CLS)
classCols[pns %in% phaseI] <- cols[1]
classCols[pns %in% phaseII.out] <- cols[2]
classCols[pns %in% phaseII] <- cols[3]
classCols[pns %in% phaseIII] <- cols[4]
classCols[pns %in% phaseIV] <- cols[5]
classCols[pns %in% phaseNA] <- cols[6]

pdf("trans.gsplom.apr11.complete.dcor1.pdf", width=60, height=60)
par(mar=c(10,5,5,5))
TRANS.R <- splom(TRANS.CLS, TRANS.DCOR, asGlyphs=T, draw.labs=T, reorder=T, grid=F, MAX=1, MIN=0.2, clust.meth="complete", main="All expressed transcription factors, rasterized, no grid, dcor weight 1", DCOR.weight=1, useRaster=T)
dev.off()

