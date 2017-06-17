#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=shared
#SBATCH --mem=5GB
#SBATCH --time=9:00:00
#SBATCH --job-name=Psb
#SBATCH --output=P%j.out

APP=$SCRATCH/fftgal/Psb
Ng=512
L=2560
Nsub=4
wisdom=${Ng}.wsdm
alpha=0.02
dK=0.01
indir=/project/projectdirs/boss/galaxy/QPM/dr12d_cubic_mocks
a=0.6452
outdir=$SCRATCH/ssm.d

echo ${SLURM_JOB_ID} starting $(date) on $(hostname)
module load gcc fftw gsl
make Psb
for id in $@
do
    log=$outdir/a${a}_$(printf '%04d' $id)/Psb.log
    time GSL_RNG_SEED=$id $APP $Ng $L $Nsub $wisdom $alpha $dK $indir $a $id $outdir 2> $log
done
echo ${SLURM_JOB_ID} ending $(date)
