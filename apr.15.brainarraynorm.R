## BrainArray annotations
## http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/16.0.0/refseq.asp
# download.file("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/16.0.0/refseq.download/celeganscerefseq.db_16.0.0.tar.gz", "celeganscerefseq.db_16.0.0.tar.gz")
# download.file("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/16.0.0/refseq.download/celeganscerefseqprobe_16.0.0.tar.gz", "celeganscerefseqprobe_16.0.0.tar.gz")
# install.packages("celeganscerefseq.db_16.0.0.tar.gz", repos=NULL, type="source")
# install.packages("celeganscerefseqprobe_16.0.0.tar.gz", repos=NULL, type="source")
library(celeganscerefseq.db)
library(celeganscerefseqprobe)
# ls("package:celeganscerefseq.db")
#Caenorhabditis_elegans
# 78.5% of 21787 probes in set

library(SCAN.UPC)


fname.ptn <- "/fs/lustre/share/yates.c.elegans.sp13/GSE2180_RAW/*.CEL.gz"
out.fpath <- "/fs/lustre/share/yates.c.elegans.sp13/"

GSE2180.SCAN = SCAN(fname.ptn, probeSummaryPackage=celeganscerefseqprobe, outFilePath=paste0(out.fpath, "GSE2180.BRAINARRAY.REFSEQ.SCAN.txt"))
#GSE2180.UPC = UPC(fname.ptn, probeSummaryPackage=celeganscerefseqprobe, outFilePath=paste0(out.fpath, "GSE2180.BRAINARRAY.REFSEQ.UPC.txt"))

save(GSE2180.SCAN, file=paste0(out.fpath,"GSE2180.BRAINARRAY.REFSEQ.SCAN.RData"))