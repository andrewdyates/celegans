## Download annotation packages
download.file("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/16.0.0/entrezg.download/celegansceentrezgprobe_16.0.0.tar.gz", "celegansceentrezgprobe_16.0.0.tar.gz")
download.file("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/16.0.0/entrezg.download/celegansceentrezg.db_16.0.0.tar.gz", "celegansceentrezg.db_16.0.0.tar.gz")

source("http://bioconductor.org/biocLite.R")
biocLite("org.Ce.eg.db")
biocLite("pd.celegans")

## optional other packages
# biocLite("celegansprobe")
# biocLite("celegans.db")
# CDF files describe which probes are part of which probe set.
# http://masker.nci.nih.gov/ev/

install.packages("celegansceentrezgprobe_16.0.0.tar.gz", repos=NULL, type="source")
install.packages("celegansceentrezg.db_16.0.0.tar.gz", repos=NULL, type="source")

library(celegansceentrezgprobe)
library(celegansceentrezg.db)

ls("package:celegansceentrezgprobe")
ls("package:celegansceentrezg.db")
