# bucki paths
library(celegansceentrezg.db)
library(energy)
library(Biobase)
source("~/code/dependency_glyph_splom/lib.R")
load("../jun5.E.GSE2180.ALL.RData")
gold.genes <- c("pal-1", "tbx-8", "tbx-9", "elt-1", "cwn-1", "hnd-1", "scrt-1", "elt-3", "nhr-25", "nob-1", "vab-7", "unc-120", "hlh-1", "mab-21", "lin-26")
load("../jun5.GSE2180.SCAN.select.tab.dcor.tab.RData")
DCOR <- M
load("../jun5.GSE2180.SCAN.select.tab.b0.0880.z2.00.bool.tab.RData")
BOOL <- M
load("../jun5.GSE2180.SCAN.select.tab.err2.th0.2000.weak.tab.RData")
WEAK <- M
load("~/celegans/GSE2180/normed/GSE2180.SCANUPC.RData")
load("~/celegans/jun5.E.GSE2180.select.RData")

IDS <- paste0(rownames(E.GSE2180.ALL.SCAN),"_at")
SYMS <- mget(IDS, celegansceentrezgSYMBOL, ifnotfound=NA)
sum(is.na(SYMS))
#[1] 0
# get gold entrez ids
gold.entrez <- rownames(E.GSE2180.ALL.SCAN)[SYMS %in% gold.genes]


# ----------------------------------------
# ALL TRANSCRIPTION FACTORS
# ----------------------------------------
trans.entrez <- read.table("~/code/c.elegans_transcription_factors/entrez_transcriptome_list.tab", sep='\t', header=T, stringsAsFactors=F)
dim(trans.entrez)
#[1] 1977    3
names(trans.entrez)
#[1] "EntrezID"            "EnsemblTranscriptID" "GeneSymbols"

stopifnot(all(rownames(DCOR)==rownames(BOOL)))
stopifnot(all(rownames(DCOR)==rownames(WEAK)))
stopifnot(all(rownames(DCOR)==colnames(DCOR)))
stopifnot(all(rownames(DCOR)==rownames(exprs(E.GSE2180.SCAN))))
stopifnot(all(rownames(DCOR) %in% rownames(E.GSE2180.ALL.SCAN)))

length(unique(trans.entrez$EntrezID)) #1254
tfentrez.plus.gold <- unique(c(trans.entrez$EntrezID, gold.entrez))
length(tfentrez.plus.gold) # 1256
qq <- rownames(DCOR) %in% tfentrez.plus.gold # select rows by entrez ID
sum(qq) # 759

DCOR.TF <- DCOR[qq,qq]
BOOL.TF <- BOOL[qq,qq]
WEAK.TF <- WEAK[qq,qq]
M.TF <- exprs(E.GSE2180.SCAN)[qq,]

write.table(BOOL.TF, file="../jun5.GSE2180.SCAN.select.tab.b0.0880.z2.00.bool.tab.TFgold.tab", sep="\t")
write.table(WEAK.TF, file="../jun5.GSE2180.SCAN.select.tab.err2.th0.2000.weak.tab.TFgold.tab", sep="\t")
write.table(DCOR.TF, file="../jun5.GSE2180.SCAN.select.tab.dcor.tab.TFgold.tab", sep="\t")
write.table(M.TF, file="../jun5.GSE2180.SCAN.select.tab.M.TFgold.tab", sep="\t")

# ------------------------------
BOOL.TF.D <- as.dist(as.matrix(read.table("../jun5.GSE2180.SCAN.select.tab.b0.0880.z2.00.bool.tab.TFgold.tab.booldist.tab", sep="\t", header=T, row.names=1, check.names=F)))
WEAK.TF.D <- as.dist(as.matrix(read.table("../jun5.GSE2180.SCAN.select.tab.err2.th0.2000.weak.tab.TFgold.tab.booldist.tab", sep="\t", header=T, row.names=1, check.names=F)))

# what is the relationship between symbols in the TF and library listings?
qq <- match(rownames(DCOR.TF), trans.entrez$EntrezID)
my.syms <- trans.entrez$GeneSymbols[qq]
brain.syms <- mget(paste0(rownames(DCOR.TF),"_at"), celegansceentrezgSYMBOL, ifnotfound=NA)
paste(my.syms, brain.syms, sep=" ::: ")
# ok, we want library symbols, not symbols in the TF mapping
stopifnot(all(!is.na(brain.syms)))
stopifnot(all(rownames(DCOR.TF)==colnames(DCOR.TF)))
stopifnot(all(rownames(BOOL.TF)==colnames(BOOL.TF)))
stopifnot(all(rownames(WEAK.TF)==colnames(WEAK.TF)))
stopifnot(all(rownames(DCOR.TF)==rownames(BOOL.TF)))
stopifnot(all(rownames(DCOR.TF)==rownames(WEAK.TF)))
# gsplom uses DCOR rownames as labels
rownames(DCOR.TF) <- brain.syms

pdf("~/GSE2180.z2.gsplom.tf.jun27.pdf", width=100, height=100)
R.GSE2180.TF <- splom(BOOL.TF, DCOR.TF, useRaster=T, sym=T, high.sat=T, MIN=0.2, MAX=1, row.cls.dist=BOOL.TF.D)
dev.off()

pdf("~/GSE2180.z2.gsplom.tf.dendro.jun27.pdf", width=180, height=12)
plot(R.GSE2180.TF$Rhclust, main="GSE2180 TF")
dev.off()

writeLines(rownames(DCOR.TF), "~/jun27.gse2180.z2.tf.gsplom.symbols.txt")
writeLines(colnames(DCOR.TF), "~/jun27.gse2180.z2.tf.gsplom.entrez.txt")
save(syms,BOOL.TF,WEAK.TF,DCOR.TF,R.GSE2180.TF,M.TF,BOOL.TF.D,WEAK.TF.D, file="../jun27.GSE2180.gsplom.RData")

pdf("~/GSE2180.z2.gsplom.tf.summary.pdf")
H.GSE2180.TF <- summary.plots(BOOL.TF, DCOR.TF, sym=T)
dev.off()

