#PBS -S /bin/bash
#PBS -q batch
#PBS -N methylC_analysis
#PBS -l nodes=1:ppn=2:HIGHMEM
#PBS -l walltime=480:00:00
#PBS -l mem=100gb

cd $PBS_O_WORKDIR
echo "Starting"
module load anaconda/3-2.2.0
module load R/3.3.1
mkdir results figures_tables
python3.4 ../../../scripts/methylC_analysis.py XXXX
Rscript ../../../scripts/methylC_analysis.R XXXX
