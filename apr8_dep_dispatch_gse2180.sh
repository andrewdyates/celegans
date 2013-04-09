#!/bin/bash
## Run dCOR, PCC, COV plus permutation tests on filtered c.elegans GSE2180

# GOLD: permute only
SRC=/nfs/01/osu6683/c.elegans/celegans.apr8.gold.tab
DST=/fs/lustre/osu6683/gse2180_gold_apr8_dep
mkdir -p $DST
cd $DST

python $HOME/pymod/dependency_matrix/permutation_test_dispatch_script.py n_permutes=3 fname=$SRC computers=[\"PCC\",\"Dcor\"] outdir=$DST n_nodes=1 n_ppn=12 hours=8
python $HOME/pymod/dependency_matrix/dispatch_script.py fname=$SRC computers=[\"PCC\",\"Dcor\"] outdir=$DST n_nodes=1 n_ppn=12 hours=1

# ------------------------------
SRC=/nfs/01/osu6683/c.elegans/celegans.apr8.expr.trans.tab
DST=/fs/lustre/osu6683/gse2180_exprs.trans_apr8_dep
mkdir -p $DST
cd $DST

python $HOME/pymod/dependency_matrix/dispatch_script.py fname=$SRC computers=[\"PCC\",\"Dcor\"] outdir=$DST n_nodes=10 n_ppn=12 hours=8
python $HOME/pymod/dependency_matrix/permutation_test_dispatch_script.py n_permutes=3 fname=$SRC computers=[\"PCC\",\"Dcor\"] outdir=$DST n_nodes=10 n_ppn=12 hours=8

# ------------------------------
SRC=/nfs/01/osu6683/c.elegans/celegans.apr8.expr.tab
DST=/fs/lustre/osu6683/gse2180_exprs_apr8_dep
mkdir -p $DST
cd $DST

python $HOME/pymod/dependency_matrix/dispatch_script.py fname=$SRC computers=[\"PCC\",\"Dcor\"] outdir=$DST n_nodes=12 n_ppn=12 hours=12
python $HOME/pymod/dependency_matrix/permutation_test_dispatch_script.py n_permutes=3 fname=$SRC computers=[\"PCC\",\"Dcor\"] outdir=$DST n_nodes=12 n_ppn=12 hours=12

