#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=shared
#SBATCH --time=4:00:00
#SBATCH --job-name=Delta
#SBATCH --output=Delta%j.out

Delta=$HOME/ssm/fftgal/Delta
Ng=512
L=2560
wisdom=${Ng}.wsdm
Nsb=4
catdir=/project/projectdirs/boss/galaxy/QPM/dr12d_cubic_mocks
a=0.6452
outdir=$SCRATCH/ssm/ana

echo ${SLURM_JOB_ID} starting $(date) on $(hostname)
module load fftw gcc
make Delta
for catid in $@
do
    log=$outdir/a${a}_$(printf '%04d' $catid)/Delta.log
    time $Delta $Ng $L $wisdom $Nsb $catdir $a $catid $outdir 2> $log
done
echo ${SLURM_JOB_ID} ending $(date)
