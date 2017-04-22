#PBS -S /bin/bash
#PBS -q batch
#PBS -N DMR_analysis
#PBS -l nodes=1:ppn=22:HIGHMEM
#PBS -l walltime=480:00:00
#PBS -l mem=100gb
#PBS -j oe

cd $PBS_O_WORKDIR
echo "Starting"
module load python/2.7.8
module load R/3.1.2
module load gsl/1.16/gcc/4.4.7

#CG
cd CG
echo "Calling CG DMRs"
python ../../../scripts/DMR_scripts/CG.py > CG_dmrfind_info.txt
cd ../
