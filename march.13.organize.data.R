source("celegans.lib.R")
load("../GSE2180.ALL.NOFILT.RData")
load("../GSE9665.ALL.NOFILT.RData")
load("../geo.GSE2180.GSE9665.RData")
# Map GSE2180 geo data to SCAN data
all(rownames(exprs(GSE2180.ALL.SCAN)) == rownames(exprs(GSE2180.ALL.UPC))) # OK
all(rownames(exprs(GSE2180.ALL.SCAN)) == rownames(exprs(GSE2180))) # FALSE
qq <- match(rownames(exprs(GSE2180)), rownames(exprs(GSE2180.ALL.SCAN)))
all(qq == 1:dim(GSE2180)[1])
any(is.na(qq)) # OK
GSE2180.ALL.SCAN <- GSE2180.ALL.SCAN[qq,]
GSE2180.ALL.UPC <- GSE2180.ALL.UPC[qq,]
all(rownames(exprs(GSE2180.ALL.SCAN)) == rownames(exprs(GSE2180))) # OK
all(featureNames(GSE2180) == rownames(exprs(GSE2180))) # OK

colnames(exprs(GSE2180.ALL.SCAN)) <- sub(".CEL.gz", "", colnames(exprs(GSE2180.ALL.SCAN)))
colnames(exprs(GSE2180.ALL.UPC)) <- sub(".CEL.gz", "", colnames(exprs(GSE2180.ALL.UPC)))

all(colnames(exprs(GSE2180.ALL.SCAN)) == colnames(exprs(GSE2180.ALL.UPC))) # OK
all(colnames(exprs(GSE2180.ALL.SCAN)) == colnames(exprs(GSE2180))) # OK
all(colnames(exprs(GSE2180.ALL.SCAN)) == sampleNames(GSE2180)) # OK

low.qual.arrays <- c("GSM39507")
col.filt <- !colnames(exprs(GSE2180)) %in% low.qual.arrays
GSE2180.GEO <- GSE2180[,col.filt]
GSE2180.SCAN <- GSE2180.ALL.SCAN[,col.filt]
GSE2180.UPC <- GSE2180.ALL.UPC[,col.filt]

featureData(GSE2180.SCAN) <- featureData(GSE2180.GEO)
featureData(GSE2180.UPC) <- featureData(GSE2180.GEO)
phenoData(GSE2180.SCAN) <- phenoData(GSE2180.GEO)
phenoData(GSE2180.UPC) <- phenoData(GSE2180.GEO)
G <- GSE2180.SCAN

## Get gold networks, add directly to annotation.
## ========================================
featureData(G)$gold<-NA

gold.row.nums <- sapply(gold.genes, function(s) grep(paste0("(",s," |",s,"$)"), featureData(GSE2180)$Gene.Symbol))
known.noexpress <- c("193640_s_at", "175709_at", "193846_at")
gold.probes <- list()
for (name in names(gold.row.nums)) {
  ids <- gold.row.nums[[name]]
  for (x in ids) {
    pname <- featureNames(G)[x]
    if (!pname %in% known.noexpress) {
      label <- paste(name, pname)
      gold.probes[label] <- x
      featureData(G)$gold[x] <- label
    }
  }
}
GSE2180.SCAN <- G; rm(G)
featureData(GSE2180.UPC)$gold <- featureData(GSE2180.SCAN)$gold
featureData(GSE2180.GEO)$gold <- featureData(GSE2180.SCAN)$gold
 
## Visually inspect, get correct probe IDs
sapply(gold.probes, function(q) featureData(GSE2180)$Gene.Symbol[q]) ## OK
sum(!is.na(featureData(GSE2180.GEO)$gold)) == 19 ## OK

save(GSE2180.SCAN, GSE2180.UPC, GSE2180.GEO, file="../GSE2180.MAR13.GOLDVAR.RData")
