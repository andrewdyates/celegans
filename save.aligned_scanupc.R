source("celegans.lib.R")
load("../GSE2180.ALL.NOFILT.RData")
load("../GSE9665.ALL.NOFILT.RData")
load("../geo.GSE2180.GSE9665.RData")
# Map GSE2180 geo data to SCAN data
all(rownames(exprs(GSE2180.ALL.SCAN)) == rownames(exprs(GSE2180.ALL.UPC))) # OK
all(rownames(exprs(GSE2180.ALL.SCAN)) == rownames(exprs(GSE2180))) # FALSE
qq <- match(rownames(exprs(GSE2180)), rownames(exprs(GSE2180.ALL.SCAN)))
all(qq == 1:dim(GSE2180)[1])
any(is.na(qq))
GSE2180.ALL.SCAN <- GSE2180.ALL.SCAN[qq,]
GSE2180.ALL.UPC <- GSE2180.ALL.UPC[qq,]
all(rownames(exprs(GSE2180.ALL.SCAN)) == rownames(exprs(GSE2180))) # OK


colnames(exprs(GSE2180.ALL.SCAN)) <- sub(".CEL.gz", "", colnames(exprs(GSE2180.ALL.SCAN)))
colnames(exprs(GSE2180.ALL.UPC)) <- sub(".CEL.gz", "", colnames(exprs(GSE2180.ALL.UPC)))

all(colnames(exprs(GSE2180.ALL.SCAN)) == colnames(exprs(GSE2180.ALL.UPC))) # OK
all(colnames(exprs(GSE2180.ALL.SCAN)) == colnames(exprs(GSE2180))) # OK

low.qual.arrays <- c("GSM39507")
col.filt <- !colnames(exprs(GSE2180)) %in% low.qual.arrays
GSE2180 <- GSE2180[,col.filt]
GSE2180.ALL.SCAN <- GSE2180.ALL.SCAN[,col.filt]
GSE2180.ALL.UPC <- GSE2180.ALL.UPC[,col.filt]

featureData(GSE2180.ALL.SCAN) <- featureData(GSE2180)
featureData(GSE2180.ALL.UPC) <- featureData(GSE2180)
phenoData(GSE2180.ALL.SCAN) <- phenoData(GSE2180)
phenoData(GSE2180.ALL.UPC) <- phenoData(GSE2180)



save(GSE2180.ALL.SCAN, GSE2180.ALL.UPC, file="../GSE2180.ALL.NOFILT.SCANUPC.CLEAN.GOLD.RData")
