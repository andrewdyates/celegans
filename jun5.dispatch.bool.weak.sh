#PBS -N boolweak2180
#PBS -l nodes=1:ppn=8
#PBS -m bea
#PBS -l walltime=4:00:00
#PBS -M yates.115.osu@gmail.com

set -x
cd $HOME
source .bash_profile

FNAME=$HOME/celegans/jun5.GSE2180.SCAN.select.tab
Z=0.27
B=0.08797455
ERR=2
time python $HOME/code/boolean_implication_fit_py/script.py fname=$FNAME b=$B z_th=$Z
time python $HOME/code/boolean_implication_fit_py/script_weak.py fname=$FNAME err=$ERR
