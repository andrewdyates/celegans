library(SCAN.UPC)
library(celegansceentrezgprobe)

fname.ptn <- "/fs/lustre/share/yates.c.elegans.sp13/GSE9965_RAW/*.cel.gz"
out.fpath <- "/fs/lustre/share/yates.c.elegans.sp13/"
GSE9665.ALL.SCAN = SCAN(fname.ptn, probeSummaryPackage=celegansceentrezgprobe, outFilePath=paste0(out.fpath, "GSE9965.ALL.SCAN.txt"))
GSE9665.ALL.UPC = UPC(fname.ptn, probeSummaryPackage=celegansceentrezgprobe, outFilePath=paste0(out.fpath, "GSE9965.ALL.UPC.txt"))
save(GSE9965.ALL.SCAN, GSE9965.ALL.UPC, file=paste0(out.fpath,"GSE9965.ALL.RData"))
