#!/bin/bash --login
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100GB
#SBATCH --job-name gene_methylation
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
#Set tmp directories
export TMPDIR=$PBS_O_WORKDIR
export TMP=$PBS_O_WORKDIR
export TEMP=$PBS_O_WORKDIR

#variables
sample=$(pwd | sed s/.*data\\/// | sed s/\\/.*//)

#get total weighted mC
echo "Gene metaplot data of $sample"
cd methylCseq
python ../../../scripts/gene_methylation.py $sample

