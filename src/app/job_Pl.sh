#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=shared
#SBATCH --mem=3GB
#SBATCH --time=3:00:00
#SBATCH --job-name=Pl
#SBATCH --output=Pl%j.out

APP=$SCRATCH/fftgal/Pl
Ng=512
L=2560
dK=0.01
wisdom=${Ng}.wsdm
fold=1
catdir=/project/projectdirs/boss/galaxy/QPM/dr12d_cubic_mocks
a=0.6452
outdir=$SCRATCH/ssm.d

echo ${SLURM_JOB_ID} starting $(date) on $(hostname)
module load gcc fftw gsl
make Pl
for catid in $@
do
    log=$outdir/a${a}_$(printf '%04d' $catid)/Pl.log
    time $APP $Ng $L $dK $wisdom $fold $catdir $a $catid $outdir 2> $log
done
echo ${SLURM_JOB_ID} ending $(date)
