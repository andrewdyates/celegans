load("../jul1.NAfiltered.tf.RData")
load("../jun27.GSE2180.gsplom.RData")
library(celegansceentrezg.db)
# for each gold gene, get a list of weak, dcor, and bool
gold.genes <- c("pal-1", "tbx-8", "tbx-9", "elt-1", "cwn-1", "hnd-1", "scrt-1", "elt-3", "nhr-25", "nob-1", "vab-7", "unc-120", "hlh-1", "mab-21", "lin-26")

#gene names saved as DCOR rownames
gold.i <- which(rownames(DCOR.TF.F) %in% gold.genes)
genenames <- rownames(DCOR.TF.F)
 ## [1] "unc-120" "cwn-1"   "hlh-1"   "mab-21"  "pal-1"   "tbx-9"   "tbx-8"  "nob-1"   "elt-1"   "nhr-25"  "elt-3"   "hnd-1"   "scrt-1"

Tabs <- list()
# note: target gene is x-axis (column)
for (i in gold.i) {
  gene <- rownames(DCOR.TF.F)[i]
  TT <- cbind(BOOL.TF.F[,i],WEAK.TF.F[,i],DCOR.TF.F[,i])
  colnames(TT) <- c("bool", "weak", "dcor")
  Tabs[[gene]] <- TT
}
# write all tables to file
for (gene in names(Tabs)) {
  write.table(Tabs[gene], file=paste0("~/Desktop/",gene,".tab"), quote=FALSE, sep="\t", col.names=NA)
}
# write tables in helen spreadsheet order (if available)


elt1.i <- readLines("~/Dropbox/biostat/local_c.elegans/entrez_orders/elt-1.txt")
elt1.names <- mget(paste0(elt1.i,"_at"), celegansceentrezgSYMBOL, ifnotfound=NA)
qq <- match(elt1.i, rownames(Tabs[["elt-1"]]))
write.table(cbind(Tabs[["elt-1"]][qq,],as.character(elt1.names)), file="~/Desktop/excel_aligned/elt-1.aligned.tab", quote=FALSE, sep="\t", col.names=NA)

hlh.i <- readLines("~/Dropbox/biostat/local_c.elegans/entrez_orders/hlh-1.txt")
hlh.names <- mget(paste0(hlh.i,"_at"), celegansceentrezgSYMBOL, ifnotfound=NA)
qq <- match(hlh.i, rownames(Tabs[["hlh-1"]]))
write.table(cbind(Tabs[["hlh-1"]][qq,],as.character(hlh.names)), file="~/Desktop/excel_aligned/hlh-1.aligned.tab", quote=FALSE, sep="\t", col.names=NA)

nob1.i <- readLines("~/Dropbox/biostat/local_c.elegans/entrez_orders/nob-1.txt")
nob1.names <- mget(paste0(nob1.i,"_at"), celegansceentrezgSYMBOL, ifnotfound=NA)
qq <- match(nob1.i, rownames(Tabs[["nob-1"]]))
write.table(cbind(Tabs[["nob-1"]][qq,],as.character(nob1.names)), file="~/Desktop/excel_aligned/nob-1.aligned.tab", quote=FALSE, sep="\t", col.names=NA)

scrt1.i <- readLines("~/Dropbox/biostat/local_c.elegans/entrez_orders/scrt-1.txt")
scrt1.names <- mget(paste0(scrt1.i,"_at"), celegansceentrezgSYMBOL, ifnotfound=NA)
qq <- match(scrt1.i, rownames(Tabs[["scrt-1"]]))
write.table(cbind(Tabs[["scrt-1"]][qq,],as.character(scrt1.names)), file="~/Desktop/excel_aligned/scrt-1.aligned.tab", quote=FALSE, sep="\t", col.names=NA)

unc120.i <- readLines("~/Dropbox/biostat/local_c.elegans/entrez_orders/unc-120.txt")
unc120.names <- mget(paste0(unc120.i,"_at"), celegansceentrezgSYMBOL, ifnotfound=NA)
qq <- match(unc120.i, rownames(Tabs[["unc-120"]]))
write.table(cbind(Tabs[["unc-120"]][qq,],as.character(unc120.names)), file="~/Desktop/excel_aligned/unc-120.aligned.tab", quote=FALSE, sep="\t", col.names=NA)

elt3.i <- readLines("~/Dropbox/biostat/local_c.elegans/entrez_orders/elt-3.txt")
elt3.names <- mget(paste0(elt3.i,"_at"), celegansceentrezgSYMBOL, ifnotfound=NA)
qq <- match(elt3.i, rownames(Tabs[["elt-3"]]))
write.table(cbind(Tabs[["elt-3"]][qq,],as.character(elt3.names)), file="~/Desktop/excel_aligned/elt-3.aligned.tab", quote=FALSE, sep="\t", col.names=NA)

hnd1.i <- readLines("~/Dropbox/biostat/local_c.elegans/entrez_orders/hnd-1.txt")
hnd1.names <- mget(paste0(hnd1.i,"_at"), celegansceentrezgSYMBOL, ifnotfound=NA)
qq <- match(hnd1.i, rownames(Tabs[["hnd-1"]]))
write.table(cbind(Tabs[["hnd-1"]][qq,],as.character(hnd1.names)), file="~/Desktop/excel_aligned/hnd-1.aligned.tab", quote=FALSE, sep="\t", col.names=NA)

nhr25.i <- readLines("~/Dropbox/biostat/local_c.elegans/entrez_orders/nhr-25.txt")
nhr25.names <- mget(paste0(nhr25.i,"_at"), celegansceentrezgSYMBOL, ifnotfound=NA)
qq <- match(nhr25.i, rownames(Tabs[["nhr-25"]]))
write.table(cbind(Tabs[["nhr-25"]][qq,],as.character(nhr25.names)), file="~/Desktop/excel_aligned/nhr-25.aligned.tab", quote=FALSE, sep="\t", col.names=NA)

pal1.i <- readLines("~/Dropbox/biostat/local_c.elegans/entrez_orders/pal-1.txt")
pal1.names <- mget(paste0(pal1.i,"_at"), celegansceentrezgSYMBOL, ifnotfound=NA)
qq <- match(pal1.i, rownames(Tabs[["pal-1"]]))
write.table(cbind(Tabs[["pal-1"]][qq,],as.character(pal1.names)), file="~/Desktop/excel_aligned/pal-1.aligned.tab", quote=FALSE, sep="\t", col.names=NA)

tbx89.i <- readLines("~/Dropbox/biostat/local_c.elegans/entrez_orders/tbx-8,9.txt")
tbx89.names <- mget(paste0(tbx89.i,"_at"), celegansceentrezgSYMBOL, ifnotfound=NA)
qq <- match(tbx89.i, rownames(Tabs[["tbx-8"]]))
write.table(cbind(Tabs[["tbx-8"]][qq,],as.character(tbx89.names)), file="~/Desktop/excel_aligned/tbx-8.aligned.tab", quote=FALSE, sep="\t", col.names=NA)
write.table(Tabs[["tbx-9"]][qq,], file="~/Desktop/excel_aligned/tbx-9.aligned.tab", quote=FALSE, sep="\t", col.names=NA)

#176450 1 3 0.7984 tbx-8
#172636 unc-120

plot(M.TF[276,],M.TF[62,], xlab="tbx-8", ylab="unc-120")
#hline(h=0.2)

which(rownames(M.TF)==176450)
#[1] 276 tbx-8
which(rownames(M.TF)==172636)
#[1] 62 unc-120

a <- M.TF[62,]>0.2  #unc-120
b <- M.TF[276,]>0.2 #tbx-8

# col unc-120 detected necessary for tbx-8 detected
sum(a&b)
#[1] 39
sum(a-(a&b)) # unc-120 detected alone
#[1] 9
sum(b-(a&b)) # tbx-8 detected alone
#[1] 2

tbx8 <- M.TF[rownames(M.TF)==176450,]
unc120 <- M.TF[rownames(M.TF)==172636,]
nob1 <- M.TF[rownames(M.TF)==176641,]
source("~/Dropbox/biostat/git_repos/boolean_implication_fit/bool.R")
source("~/Dropbox/biostat/git_repos/boolean_implication_fit/step.up.R")
R.tbx8 <- fit.upstep(tbx8)     # R.tbx8$th
R.unc120 <- fit.upstep(unc120) # R.unc120$th
R.nob1 <- fit.upstep(nob1)     # R.nob1$th
b<-0.0880

plot(tbx8,unc120)
pdf("~/Desktop/tbx8.unc120.pdf")
cls.pair(x=tbx8, y=unc120, x.th=R.tbx8$th, y.th=R.unc120$th, xlab="tbx8", ylab="unc120", Z=2, b.x=b, b.y=b, do.plot=T)
dev.off()

pdf("~/Desktop/unc120.tbx8.pdf")
cls.pair(x=unc120, y=tbx8, x.th=R.unc120$th, y.th=R.tbx8$th, xlab="unc120", ylab="tbx8", Z=2, b.x=b, b.y=b, do.plot=T)
abline(v=0.2, col="green", lwd=1)
abline(h=0.2, col="green", lwd=2)
dev.off()

pdf("~/Desktop/nob1.tbx8.pdf")
cls.pair(x=nob1, y=tbx8, x.th=R.nob1$th, y.th=R.tbx8$th, xlab="nob1", ylab="tbx8", Z=2, b.x=b, b.y=b, do.plot=T)
abline(v=0.2, col="green", lwd=1)
abline(h=0.2, col="green", lwd=2)
dev.off()


ee <- c(172844,172981,174341,174614,176020,177016,178465,178919,179399,179589,180324,180431,181229,181302,185718,187182,187186,191703,191719,172432,172728,176450,172626,172981,174341,174614,174721,175612,176241,176632,177016,177618,178919,179276,179399,180324,180848,181229,181529,181551,184557,185593,185718,186627,191149,191703,191719,172432,184793,176450)

nn <- mget(paste0(as.character(ee),"_at"), celegansceentrezgSYMBOL, ifnotfound=NA)
