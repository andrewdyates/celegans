SRC=/fs/lustre/osu6683/gse2180_exprs_apr8_dep
PKL_FNAME_PCC=$SRC/compiled_dep_matrices/celegans.apr8.expr.pkl.PEARSON.values.pkl
PKL_FNAME_DCOR=$SRC/compiled_dep_matrices/celegans.apr8.expr.pkl.DCOR.values.pkl
ROW_FNAME=$SRC/celegans.apr8.expr.rowIDs.txt
COL_FNAME=$SRC/celegans.apr8.expr.rowIDs.txt
OUTDIR=$HOME/c.elegans

echo "PCC..."
/usr/bin/time python $HOME/pymod/pkl_txt_RData/pkl_to_RData.py pkl_fname=$PKL_FNAME_PCC row_fname=$ROW_FNAME col_fname=$COL_FNAME outdir=$OUTDIR

echo "DCOR..."
/usr/bin/time python $HOME/pymod/pkl_txt_RData/pkl_to_RData.py pkl_fname=$PKL_FNAME_DCOR row_fname=$ROW_FNAME col_fname=$COL_FNAME outdir=$OUTDIR 
