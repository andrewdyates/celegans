library(SCAN.UPC)
library(celegansceentrezgprobe)

fname.ptn <- "/fs/lustre/share/yates.c.elegans.sp13/GSE2180_RAW/*.CEL.gz"
out.fpath <- "/fs/lustre/share/yates.c.elegans.sp13/"
GSE2180.ALL.SCAN = SCAN(fname.ptn, probeSummaryPackage=celegansceentrezgprobe, outFilePath=paste0(out.fpath, "GSE2180.ALL.SCAN.txt"))
GSE2180.ALL.UPC = UPC(fname.ptn, probeSummaryPackage=celegansceentrezgprobe, outFilePath=paste0(out.fpath, "GSE2180.ALL.UPC.txt"))
save(GSE2180.ALL.SCAN, GSE2180.ALL.UPC, file="fs/lustre/share/yates.c.elegans.sp13/GSE2180.ALL.RData")
