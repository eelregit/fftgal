#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=shared
#SBATCH --time=3:00:00
#SBATCH --job-name=Pfd
#SBATCH --output=P%j.out

APP=$SCRATCH/fftgal/Pl
Ng=128
L=2560
fold=4
wisdom=${Ng}.wsdm
dK=0.01
indir=/project/projectdirs/boss/galaxy/QPM/dr12d_cubic_mocks
a=0.6452
outdir=$SCRATCH/ssm.d

echo ${SLURM_JOB_ID} starting $(date) on $(hostname)
module load gcc fftw gsl
make Pl
for id in $@
do
    log=$outdir/a${a}_$(printf '%04d' $id)/Pfd.log
    time $APP $Ng $L $fold $wisdom $dK $indir $a $id $outdir 2> $log
done
echo ${SLURM_JOB_ID} ending $(date)
