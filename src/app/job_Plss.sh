#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=shared
#SBATCH --time=2:00:00
#SBATCH --job-name=Plss
#SBATCH --output=Pl%j.out

Plss=$SCRATCH/fftgal/Plss
Ng=256
L=2560
wisdom=${Ng}.wsdm
Nss=4
catdir=/project/projectdirs/boss/galaxy/QPM/dr12d_cubic_mocks
a=0.6452
outdir=$SCRATCH/ssm.d

echo ${SLURM_JOB_ID} starting $(date) on $(hostname)
module load fftw gcc
make Plss
for catid in $@
do
    log=$outdir/a${a}_$(printf '%04d' $catid)/Plss.log
    time $Plss $Ng $L $wisdom $Nss $catdir $a $catid $outdir 2> $log
done
echo ${SLURM_JOB_ID} ending $(date)
