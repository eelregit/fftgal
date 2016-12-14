#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=shared
#SBATCH --time=4:00:00
#SBATCH --job-name=Plsb
#SBATCH --output=Pl%j.out

Plsb=$HOME/ssm/fftgal/Plsb
Ng=128
L=2560
wisdom=${Ng}.wsdm
Nsb=4
catdir=/project/projectdirs/boss/galaxy/QPM/dr12d_cubic_mocks
a=0.6452
outdir=$HOME/ssm/ana

echo ${SLURM_JOB_ID} starting $(date) on $(hostname)
#module load fftw gcc
#make Plsb
for catid in $@
do
    log=$outdir/a${a}_$(printf '%04d' $catid)/Pl.log
    time $Plsb $Ng $L $wisdom $Nsb $catdir $a $catid $outdir 2> $log
done
echo ${SLURM_JOB_ID} ending $(date)
