library(SCAN.UPC)
#library(celegansceentrezgprobe)

fname.ptn <- "/fs/lustre/share/yates.c.elegans.sp13/GSE9665_RAW/*.cel.gz"
out.fpath <- "/fs/lustre/share/yates.c.elegans.sp13/"

GSE9665.ALL.SCAN = SCAN(fname.ptn, outFilePath=paste0(out.fpath, "GSE9665.ALL.SCAN.NOFILT.txt"))
GSE9665.ALL.UPC = UPC(fname.ptn, outFilePath=paste0(out.fpath, "GSE9665.ALL.UPC.NOFILT.txt"))
save(GSE9665.ALL.SCAN, GSE9665.ALL.UPC, file=paste0(out.fpath,"GSE9665.ALL.NOFILT.RData"))
