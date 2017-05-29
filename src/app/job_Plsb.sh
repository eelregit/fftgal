#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=shared
#SBATCH --time=2:00:00
#SBATCH --job-name=Plsb
#SBATCH --output=Pl%j.out

Plsb=$SCRATCH/fftgal/Plsb
Ng=128
L=2560
wisdom=${Ng}.wsdm
Nsb=4
catdir=/project/projectdirs/boss/galaxy/QPM/dr12d_cubic_mocks
a=0.6452
outdir=$SCRATCH/ssm.d

echo ${SLURM_JOB_ID} starting $(date) on $(hostname)
module load gcc fftw gsl
make Plsb
for catid in $@
do
    log=$outdir/a${a}_$(printf '%04d' $catid)/Plsb.log
    time $Plsb $Ng $L $wisdom $Nsb $catdir $a $catid $outdir 2> $log
done
echo ${SLURM_JOB_ID} ending $(date)
