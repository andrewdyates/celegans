library(SCAN.UPC)
#library(celegansceentrezgprobe)

fname.ptn <- "/fs/lustre/share/yates.c.elegans.sp13/GSE2180_RAW/GSM39425.CEL.gz"
out.fpath <- "/fs/lustre/share/yates.c.elegans.sp13/"
M1 = SCAN(fname.ptn, outFilePath=paste0(out.fpath, "GSE2180.TEST1.SCAN.txt"))

library(celegansceentrezgprobe)
M2 = SCAN(fname.ptn, probeSummaryPackage=celegansceentrezgprobe, outFilePath=paste0(out.fpath, "GSE2180.TEST2.SCAN.txt"))