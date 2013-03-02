wget http://www.ncbi.nlm.nih.gov/geosuppl/?acc=GSE2180
wget http://www.ncbi.nlm.nih.gov/geosuppl/?acc=GSE9665
OUTDIR=/fs/lustre/share/yates.c.elegans.sp13
python $HOME/pymod/geo_downloader/script.py gse_id=GSE2180 outdir=$OUTDIR
python $HOME/pymod/geo_downloader/script.py gse_id=GSE9665 outdir=$OUTDIR

