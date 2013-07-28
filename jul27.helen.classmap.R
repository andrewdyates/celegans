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
hlh.names <- mget(paste0(elt1.i,"_at"), celegansceentrezgSYMBOL, ifnotfound=NA)
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

