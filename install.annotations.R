download.file("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/16.0.0/entrezg.download/celegansceentrezgprobe_16.0.0.tar.gz", "celegansceentrezgprobe_16.0.0.tar.gz")
download.file("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/16.0.0/entrezg.download/celegansceentrezg.db_16.0.0.tar.gz", "celegansceentrezg.db_16.0.0.tar.gz")

source("http://bioconductor.org/biocLite.R")
biocLite("org.Ce.eg.db")

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


## > install.packages("celegansceentrezg.db_16.0.0.tar.gz", repos=NULL, type="source")
## * installing *source* package ‘celegansceentrezg.db’ ...
## ** R
## ** inst
## ** preparing package for lazy loading
## ** help
## *** installing help indices
## ** building package indices
## ** testing if installed package can be loaded
## Error : .onLoad failed in loadNamespace() for 'celegansceentrezg.db', details:
##   call: sqliteExecStatement(con, statement, bind.data)
##   error: RS-DBI driver: (error in statement: no such table: chrlengths)
## Error: loading failed
## Execution halted
## ERROR: loading failed
## * removing ‘/nfs/01/osu6683/R/R-2.15.1/library/celegansceentrezg.db’
## Warning message:
## In install.packages("celegansceentrezg.db_16.0.0.tar.gz", repos = NULL,  :
##   installation of package ‘celegansceentrezg.db_16.0.0.tar.gz’ had non-zero exit status