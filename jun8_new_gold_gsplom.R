# bucki paths
library(celegansceentrezg.db)
source("~/code/dependency_glyph_splom/lib.R")
gold.genes <- c("pal-1", "tbx-8", "tbx-9", "elt-1", "cwn-1", "hnd-1", "scrt-1", "elt-3", "nhr-25", "nob-1", "vab-7", "unc-120", "hlh-1", "mab-21", "lin-26")
load("../jun5.GSE2180.SCAN.select.tab.dcor.tab.RData")
DCOR <- M
load("../jun5.GSE2180.SCAN.select.tab.b0.0880.z0.27.bool.tab.RData")
BOOL <- M
load("../jun5.GSE2180.SCAN.select.tab.err2.th0.2000.weak.tab.RData")
WEAK <- M
load("~/celegans/GSE2180/normed/GSE2180.SCANUPC.RData")
load("~/celegans/jun5.E.GSE2180.select.RData")


E.GSE2180.SCAN
IDS <- paste0(rownames(DCOR),"_at")
SYMS <- mget(IDS, celegansceentrezgSYMBOL, ifnotfound=NA)
qq <- na.omit(match(gold.genes, SYMS))


GOLD.M <- exprs(E.GSE2180.SCAN)[qq,]
rownames(GOLD.M) <- SYMS[qq]
write.table(GOLD.M, file="~/new.gold.celegans.tab", sep="\t", quote=F)

SYMS[qq[1]]
SYMS[qq[2]]
pal1 <- exprs(E.GSE2180.SCAN)[qq[1],]
tbx8 <- exprs(E.GSE2180.SCAN)[qq[2],]
dcor(exprs(E.GSE2180.SCAN)[qq[1],],exprs(E.GSE2180.SCAN)[qq[2],])
dcor(exprs(E.GSE2180.SCAN)[qq[1],],exprs(E.GSE2180.SCAN)[qq[3],])

DCOR.gold <- DCOR[qq,qq]
BOOL.gold <- BOOL[qq,qq]
WEAK.gold <- WEAK[qq,qq]
rownames(DCOR.gold) <- SYMS[qq]
colnames(DCOR.gold) <- SYMS[qq]
rownames(BOOL.gold) <- SYMS[qq]
colnames(BOOL.gold) <- SYMS[qq]

pdf("~/jun8.newgold.gsplom.pdf")
R.GOLD <- splom(BOOL.gold, DCOR.gold, useRaster=F, MIN=0.32, MAX=1, grid=F, reorder=T)
dev.off()
#SYMS.ALL <- mget(paste0(rownames(S.GSE2180),"_at"), celegansceentrezgSYMBOL, ifnotfound=NA)
#match(gold.genes, SYMS.ALL) # lin-26 (3565051) is missing

# TODO: all transcription factors
