#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=shared
#SBATCH --time=3:00:00
#SBATCH --job-name=DLsb
#SBATCH --output=DL%j.out

APP=$SCRATCH/fftgal/DLsb
Ng=512
L=2560
wisdom=${Ng}.wsdm
Nsub=4
catdir=/project/projectdirs/boss/galaxy/QPM/dr12d_cubic_mocks
a=0.6452
outdir=$SCRATCH/ssm.d

echo ${SLURM_JOB_ID} starting $(date) on $(hostname)
module load gcc fftw gsl
make DLsb
for catid in $@
do
    log=$outdir/a${a}_$(printf '%04d' $catid)/DLsb.log
    time $APP $Ng $L $wisdom $Nsub $catdir $a $catid $outdir 2> $log
done
echo ${SLURM_JOB_ID} ending $(date)
