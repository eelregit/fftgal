#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=shared
#SBATCH --time=4:00:00
#SBATCH --job-name=Pl
#SBATCH --output=Pl%j.out

Pl=$SCRATCH/fftgal/Pl
Ng=512
L=2560
wisdom=${Ng}.wsdm
catdir=/project/projectdirs/boss/galaxy/QPM/dr12d_cubic_mocks
a=0.6452
outdir=$SCRATCH/ssm.d

echo ${SLURM_JOB_ID} starting $(date) on $(hostname)
module load fftw gcc
make Pl
for catid in $@
do
    log=$outdir/a${a}_$(printf '%04d' $catid)/Pl.log
    time $Pl $Ng $L $wisdom $catdir $a $catid $outdir 2> $log
done
echo ${SLURM_JOB_ID} ending $(date)
