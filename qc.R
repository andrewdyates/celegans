library(sva)
library(lumi)
library(celegansceentrezg.db)

load("../GSE2180.ALL.RData")
load("../GSE9665.ALL.RData")

## Load GSE2180 phenotype data.
## ------------------------------
AttrT <- read.table("../GSE2180_GPL200.samples.tab", sep="\t", header=TRUE, row.names=1)
AttrT <- as.data.frame(t(AttrT))
AttrT$time <- as.numeric(as.character(AttrT$"n:Sample_time"))
AttrT$genotype <- AttrT$"n:Sample_genotype"
# Title information in AttrT$title?
sampleNames(GSE2180.ALL.SCAN) <- sub(".CEL.gz", "", sampleNames(GSE2180.ALL.SCAN))
sampleNames(GSE2180.ALL.UPC) <- sub(".CEL.gz", "", sampleNames(GSE2180.ALL.UPC))
pData(phenoData(GSE2180.ALL.SCAN)) <- AttrT
pData(phenoData(GSE2180.ALL.UPC)) <- AttrT


## Load GSE9665 phenotype data.
## ------------------------------
AttrT <- read.table("../GSE9665_GPL200.samples.tab", sep="\t", header=TRUE, row.names=1)
AttrT <- as.data.frame(t(AttrT))
AttrT$knockout <- as.factor(sub("mex-3 (\\w+-[0-9,]+)? ?(\\d+) (\\w)", AttrT$title, replacement="\\1"))
AttrT$num <- as.factor(sub("mex-3 (\\w+-[0-9,]+)? ?(\\d+) (\\w)", AttrT$title, replacement="\\2"))
AttrT$primer <- as.factor(sub("mex-3 (\\w+-[0-9,]+)? ?(\\d+) (\\w).*", AttrT$title, replacement="\\3"))

sampleNames(GSE9665.ALL.SCAN) <- sub(".cel.gz", "", sampleNames(GSE9665.ALL.SCAN))
sampleNames(GSE9665.ALL.UPC) <- sub(".cel.gz", "", sampleNames(GSE9665.ALL.UPC))
pData(phenoData(GSE9665.ALL.SCAN)) <- AttrT
pData(phenoData(GSE9665.ALL.UPC)) <- AttrT

## > dim(GSE2180.ALL.SCAN)
## Features  Samples 
##    17079      123 
## > dim(GSE9665.ALL.SCAN)
## Features  Samples 
##    17079       74

inspect <- function(M.SCAN, M.UPC, name="GSE2180", th.probe=0.6, th.array=0.21) {
  upc.max1 <- apply(exprs(M.UPC), 1, max)
  upc.filt.med1 <- apply(exprs(M.UPC[upc.max1>th.probe,]), 1, median)
  row.filt <- upc.max1>th.probe
  print(length(row.filt)); print(sum(row.filt))
  
  pdf(paste0(name,".hist.probe.UPC.pdf"), width=8, height=8)
  hist(upc.max1, main=paste0(name," All probe max UPC"))
  hist(upc.filt.med1, main=paste0(name," Median UPC Filtered probe (max >",th.probe,")"))
  dev.off()

  upc.filt.means2 <- apply(exprs(M.UPC[row.filt,]), 2, mean)
  col.filt <- upc.filt.means2>=th.array
  print(length(col.filt)); print(sum(col.filt)); 
  
  print(sampleNames(M.UPC)[!col.filt])
  ## [1] "GSM39482" "GSM39505" "GSM39507"

  n<-length(upc.max1)
  col<-rep('white',n)
  col[!col.filt]<-'red'

  pdf(paste0(name,".boxplots.probefilt.pdf"), width=20, height=8)
  boxplot(M.UPC[row.filt,], col=col, main=paste0(name," UPC"))
  boxplot(M.SCAN[row.filt,], col=col, main=paste0(name," SCAN"))
  dev.off()

  pdf(paste0(name,".density.probefilt.pdf"), width=8, height=12)
  density(M.UPC[row.filt,col.filt], main=paste0(name," UPC Array Good"))
  if(sum(!col.filt)>0) {
    density(M.UPC[row.filt,!col.filt], main=paste0(name," UPC Array Bad (mean<",th.array,")"))
  }
  dev.off()
  R = list()
  R$row.filt <- row.filt
  R$col.filt <- col.filt
  R
}

GSE2180.R <- inspect(GSE2180.ALL.SCAN, GSE2180.ALL.UPC, "GSE2180", 0.6, 0.21)
GSE9665.R <- inspect(GSE9665.ALL.SCAN, GSE9665.ALL.UPC, "GSE9665", 0.6, 0.21)

GSE2180 <- GSE2180.ALL.SCAN[GSE2180.R$row.filt, GSE2180.R$col.filt]
GSE9665 <- GSE9665.ALL.SCAN[GSE9665.R$row.filt, GSE9665.R$col.filt]

test.var.diff <- function(var, M) {
  # Return fraction of rows with significant diff expression with batch variable.
  mod <- model.matrix(~var)
  mod0 <- cbind(mod[,1])
  pp <- p.adjust(f.pvalue(M,mod,mod0), method="BH")
  mean(pp < 0.05)
}
test.var.diff(GSE9665$primer, exprs(GSE9665))
[1] 0


## EXPORT DATA
save(GSE2180, file="../GSE2180.clean.RData")
save(GSE9665, file="../GSE9665.clean.RData")
write.table(exprs(GSE2180), file="../GSE2180.clean.tab", quote=FALSE, sep="\t")
write.table(exprs(GSE9665), file="../GSE9665.clean.tab", quote=FALSE, sep="\t")
