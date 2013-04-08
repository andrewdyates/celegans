library(lumi) # for nice figures, plus Biobase
library(Biobase)
load("../GSE2180.MAR13.GOLDVAR.RData")
gold.i <- which(!is.na(featureData(GSE2180.SCAN)$gold))
qq.c <- GSE2180.GEO$genotype %in% c('N2','ms')

E <- GSE2180.SCAN[,qq.c]
P <- GSE2180.UPC[,qq.c]
row.maxs <- apply(exprs(E),1,max)
num.over.5 <- apply(exprs(E) > 0.5,1,sum)
upc.num.over.5 <- apply(exprs(P) > 0.5,1,sum)

png("GSE2180.SCAN.vs.UPC.N2ms.png")
plot(exprs(E), exprs(P), main="SCAN expression vs UPC percentage")
dev.off()

sum(upc.num.over.5>=4)
#[1] 7063
sum(num.over.5>=4)
#[1] 7698
upc.num.over.5[gold.i]
##   174037_at   174043_at   175663_at   175771_at   175801_at 188239_s_at 
##          55          60          17          16           4           5 
## 188486_s_at   188706_at   190309_at   190477_at 190539_s_at 191468_s_at 
##          34          15          20          10          52           8 
##   192307_at 192655_s_at   192707_at 193001_s_at   193341_at 193342_s_at 
##           7          28          28           0          37          53 
##   193759_at 
##          11 

pdf("GSE2180.num.exprs.pdf")
hist(upc.num.over.5[upc.num.over.5>=4], breaks=77, xlim=c(4,61), main="count over UPC>0.5")
dev.off()

# Filter by minimum expression confidence
expr.probes.qq <- upc.num.over.5>=4

trans.factors <- unique(as.vector(unlist(read.table("celegans.transcription.factors.txt", stringsAsFactors=F))))
trans.row.nums <- sapply(trans.factors, function(s) grep(paste0("(",s," |",s,"$)"), featureData(E)$Gene.Symbol))

# Assign gene label to expression matrix. Warning: this is inefficient.
featureData(E)$transfact.label<-NA
featureData(E)$transfact.gene<-NA
transfact.probes <- list()
for (name in names(trans.row.nums)) {
  ids <- trans.row.nums[[name]]
  for (x in ids) {
    pname <- featureNames(E)[x]
    old.label <- featureData(E)$transfact.gene[x]
    if (!is.na(old.label))
      label <- paste0(old.label, ";", name, " ", pname)
    else
      label <- paste0(name, " ", pname)
    transfact.probes[label] <- x
    featureData(E)$transfact.label[x] <- label # unique, builds on past designations
    featureData(E)$transfact.gene[x] <- name # may clobber existing gene designations
  }
}

# of 22,625 probes, 792 are associated with at least one transcription factor
sum(!is.na(featureData(E)$transfact.gene))
# [1] 792
length(unique(featureData(E)$transfact.gene))
# [1] accounting for 574 unique genes from 606 total factors




save(expr.probes.qq, M, file="probe.filtered.celegans.RData")


