#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=shared
#SBATCH --time=4:00:00
#SBATCH --job-name=ssm
#SBATCH --output=ssm%j.out

ssm=$HOME/ssm/fftgal/ssm
Ng=512
L=2560
wisdom=${Ng}.wsdm
Nsb=4
catdir=/project/projectdirs/boss/galaxy/QPM/dr12d_cubic_mocks
a=0.6452
outdir=$SCRATCH/ssm.d

echo ${SLURM_JOB_ID} starting $(date) on $(hostname)
module load fftw gcc
make ssm
for catid in $@
do
    log=$outdir/a${a}_$(printf '%04d' $catid)/ssm.log
    time $ssm $Ng $L $wisdom $Nsb $catdir $a $catid $outdir 2> $log
done
echo ${SLURM_JOB_ID} ending $(date)
