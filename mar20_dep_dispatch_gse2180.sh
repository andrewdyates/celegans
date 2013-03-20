#!/bin/bash
## Run dCOR, PCC, COV plus permutation tests on filtered c.elegans GSE2180

SRC=/nfs/01/osu6683/c.elegans/GSE2180.SCAN.N2.ms.max_0.5.tab
DST=/fs/lustre/osu6683/gse2180_mar20_dep
mkdir -p $DST
cd $DST

python $HOME/pymod/dependency_matrix/dispatch_script.py fname=$SRC computers=[\"PCC\",\"Cov\",\"Dcor\"] outdir=$DST n_nodes=10 n_ppn=12 hours=30
python $HOME/pymod/dependency_matrix/permutation_test_dispatch_script.py n_permutes=1 fname=$SRC computers=[\"PCC\",\"Dcor\"] outdir=$DST n_nodes=10 n_ppn=12 hours=30

