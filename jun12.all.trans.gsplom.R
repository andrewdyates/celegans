# bucki paths
library(celegansceentrezg.db)
library(energy)
source("~/code/dependency_glyph_splom/lib.R")
load("../jun5.E.GSE2180.ALL.RData")
gold.genes <- c("pal-1", "tbx-8", "tbx-9", "elt-1", "cwn-1", "hnd-1", "scrt-1", "elt-3", "nhr-25", "nob-1", "vab-7", "unc-120", "hlh-1", "mab-21", "lin-26")
load("../jun5.GSE2180.SCAN.select.tab.dcor.tab.RData")
DCOR <- M
load("../jun5.GSE2180.SCAN.select.tab.b0.0880.z0.27.bool.tab.RData")
BOOL <- M
load("../jun5.GSE2180.SCAN.select.tab.err2.th0.2000.weak.tab.RData")
WEAK <- M
load("~/celegans/GSE2180/normed/GSE2180.SCANUPC.RData")
load("~/celegans/jun5.E.GSE2180.select.RData")

IDS <- paste0(rownames(E.GSE2180.ALL.SCAN),"_at")
SYMS <- mget(IDS, celegansceentrezgSYMBOL, ifnotfound=NA)
# get gold entrez ids
gold.entrez <- rownames(E.GSE2180.ALL.SCAN)[SYMS %in% gold.genes]


# ----------------------------------------
# ALL TRANSCRIPTION FACTORS
# ----------------------------------------
trans.entrez <- read.table("~/code/c.elegans_transcription_factors/entrez_transcriptome_list.tab", sep='\t', header=T, stringsAsFactors=F)
dim(trans.entrez)
#[1] 1977    3

stopifnot(all(rownames(DCOR)==rownames(BOOL)))
stopifnot(all(rownames(DCOR)==rownames(WEAK)))
stopifnot(all(rownames(DCOR)==colnames(DCOR)))

length(unique(trans.entrez$EntrezID)) #1254
tfentrez.plus.gold <- unique(c(trans.entrez$EntrezID, gold.entrez)) # 1256
qq <- rownames(DCOR) %in% tfentrez.plus.gold
sum(qq) # 759

DCOR.TF <- DCOR[qq,qq]
BOOL.TF <- BOOL[qq,qq]
WEAK.TF <- WEAK[qq,qq]

write.table(BOOL.TF, file="../jun5.GSE2180.SCAN.select.tab.b0.0880.z0.27.bool.tab.TFgold.tab", sep="\t")
write.table(WEAK.TF, file="../jun5.GSE2180.SCAN.select.tab.err2.th0.2000.weak.tab.TFgold.tab", sep="\t")
write.table(DCOR.TF, file="../jun5.GSE2180.SCAN.select.tab.dcor.tab.TFgold.tab", sep="\t")

BOOL.TF.D <- as.dist(as.matrix(read.table("../jun5.GSE2180.SCAN.select.tab.b0.0880.z0.27.bool.tab.TFgold.tab.booldist.tab", sep="\t", header=T, row.names=1, check.names=F)))
WEAK.TF.D <- as.dist(as.matrix(read.table("../jun5.GSE2180.SCAN.select.tab.err2.th0.2000.weak.tab.TFgold.tab.booldist.tab", sep="\t", header=T, row.names=1, check.names=F)))

qq <- match(rownames(DCOR.TF), trans.entrez$EntrezID)
syms <- trans.entrez$GeneSymbols[qq]
qq <- grep("\\|",syms)
syms[qq] <- trans.entrez$EnsemblTranscriptID[qq]

rownames(DCOR.TF) <- syms

pdf("~/GSE2180.gsplom.tf.pdf", width=100, height=100)
R.GSE2180.TF <- splom(BOOL.TF, DCOR.TF, useRaster=T, sym=T, high.sat=T, MIN=0.2, MAX=1, row.cls.dist=BOOL.TF.D)
dev.off()

pdf("~/GSE2180.gsplom.tf.dendro.pdf", width=180, height=12)
plot(R.GSE2180.TF$Rhclust, main="GSE2180 TF")
dev.off()

save(syms,BOOL.TF,WEAK.TF,DCOR.TF,R.GSE2180.TF, file="../jun12.GSE2180.gsplom.RData")
