python ~/code/qsub/script.py script="time python ~/code/all_bool_dist/script.py fname=~/celegans/jun5.GSE2180.SCAN.select.tab.b0.0880.z0.27.bool.tab.TFgold.tab" options="#PBS -M yates.115.osu@gmail.com" jobname="31448.TF.BOOL" n_ppn=8 email=True walltime="6:00:00"

python ~/code/qsub/script.py script="time python ~/code/all_bool_dist/script.py fname=~/celegans/jun5.GSE2180.SCAN.select.tab.err2.th0.2000.weak.tab.TFgold.tab" options="#PBS -M yates.115.osu@gmail.com" jobname="2180.TF.BOOL" n_ppn=8 email=True walltime="6:00:00"


python ~/code/qsub/script.py script="time python ~/code/all_bool_dist/script.py fname=~/celegans/jun5.GSE2180.SCAN.select.tab.b0.0880.z2.00.bool.tab.TFgold.tab" options="#PBS -M yates.115.osu@gmail.com" jobname="2180.TF.BOOL" n_ppn=8 email=True walltime="6:00:00"
