library("Biobase")
load("../celegans.apr8.expr.RData")
load("../apr9.dcor.cls.RData")
source("~/Dropbox/biostat/git_repos/dependency_glyph_splom/lib.R") # local ref
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

