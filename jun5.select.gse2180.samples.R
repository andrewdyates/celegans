# bucki paths
library(Biobase)
library(celegansceentrezg.db)
load("~/celegans/GSE2180/normed/GSE2180.SCANUPC.RData")
gold.genes <- c("pal-1", "tbx-8", "tbx-9", "elt-1", "cwn-1", "hnd-1", "scrt-1", "elt-3", "nhr-25", "nob-1", "vab-7", "unc-120", "hlh-1", "mab-21", "lin-26")
# [1] "S.GSE2180" "U.GSE2180"

AttrT <- read.table("~/celegans/GSE2180/raw/GSE2180_GPL200.samples.tab", sep="\t", header=TRUE, row.names=1, stringsAsFactors=F)
AttrT <- as.data.frame(t(AttrT))
AttrT$time <- as.numeric(as.character(AttrT$"n:Sample_time"))
AttrT$genotype <- AttrT$"n:Sample_genotype"

all(rownames(S.GSE2180) == rownames(U.GSE2180)) # T
colnames(S.GSE2180) <- sub(".CEL.gz", "", colnames(S.GSE2180))
colnames(U.GSE2180) <- sub(".CEL.gz", "", colnames(U.GSE2180))
all(colnames(S.GSE2180) == colnames(U.GSE2180)) # F
# align column names to each other
qq <- match(colnames(S.GSE2180), colnames(U.GSE2180))
U.GSE2180 <- U.GSE2180[,qq]
all(colnames(S.GSE2180) == colnames(U.GSE2180)) # T

# create expression sets
E.GSE2180.ALL.SCAN <- ExpressionSet(S.GSE2180)
E.GSE2180.ALL.UPC  <- ExpressionSet(U.GSE2180)
all(sampleNames(E.GSE2180.ALL.SCAN)==sampleNames(E.GSE2180.ALL.UPC))
all(sampleNames(E.GSE2180.ALL.SCAN)==colnames(S.GSE2180))
all(sampleNames(E.GSE2180.ALL.SCAN)==colnames(exprs(E.GSE2180.ALL.SCAN)))
all(sampleNames(E.GSE2180.ALL.UPC)==colnames(U.GSE2180))
all(sampleNames(E.GSE2180.ALL.UPC)==colnames(exprs(E.GSE2180.ALL.UPC))) # OK

# align attribute table, add to exprSets
all(rownames(AttrT) == AttrT$geo_accession)
all(rownames(AttrT) == sampleNames(E.GSE2180.ALL.SCAN))
qq <- match(sampleNames(E.GSE2180.ALL.SCAN), rownames(AttrT))
AttrT <- AttrT[qq,]
all(rownames(AttrT) == sampleNames(E.GSE2180.ALL.SCAN)) # T
all(rownames(AttrT) == AttrT$geo_accession)
pData(phenoData(E.GSE2180.ALL.SCAN)) <- AttrT
pData(phenoData(E.GSE2180.ALL.UPC)) <- AttrT

# add gene symbol from brainarray to featureData
IDS <- paste(rownames(exprs(E.GSE2180.ALL.SCAN)),"at",sep="_")
SYMS <- mget(IDS, celegansceentrezgSYMBOL, ifnotfound=NA)
featureData(E.GSE2180.ALL.SCAN)$SYM <- SYMS
featureData(E.GSE2180.ALL.UPC)$SYM <- SYMS

# save aligned expr sets
save(E.GSE2180.ALL.SCAN, E.GSE2180.ALL.UPC, file="../jun5.E.GSE2180.ALL.RData")

# select wildtype and mex-3 mutants only
col.qq <- E.GSE2180.ALL.SCAN$genotype %in% c("N2", "ms")
E.GSE2180.SCAN <- E.GSE2180.ALL.SCAN[,col.qq]
E.GSE2180.UPC <- E.GSE2180.ALL.UPC[,col.qq]

# quality control: remove putatively low quality array determined by visual inspection
low.qual.arrays <- c("GSM39507") # by manual inspection
qq <- !sampleNames(E.GSE2180.SCAN) %in% low.qual.arrays
E.GSE2180.SCAN <- E.GSE2180.SCAN[,qq]
E.GSE2180.UPC <- E.GSE2180.UPC[,qq]

# get std distribution, select noise parameter
stds <- apply(exprs(E.GSE2180.SCAN),1,sd)
b <- as.numeric(quantile(stds, 0.03)*2)
# 0.08797455
quantile(stds,0.25)
# 0.06306032
# filter probesets with std less than b
# 16811 to 8770 features
E.GSE2180.SCAN <- E.GSE2180.SCAN[stds>b,]
E.GSE2180.UPC <- E.GSE2180.UPC[stds>b,]

# manual checks for status of gold genes
which(rownames(exprs(E.GSE2180.SCAN))=="175607") # mab-21: 3445
gold.qq <- which(featureData(E.GSE2180.SCAN)$SYM %in% gold.genes)
apply(exprs(E.GSE2180.SCAN[gold.qq,]),1,max)
stds[gold.qq]

# save filtered expression matrices
save(E.GSE2180.SCAN, E.GSE2180.UPC, file="../jun5.E.GSE2180.select.RData")
write.table(exprs(E.GSE2180.SCAN), file="../jun5.GSE2180.SCAN.select.tab", sep="\t", row.names=T, col.names=T)
