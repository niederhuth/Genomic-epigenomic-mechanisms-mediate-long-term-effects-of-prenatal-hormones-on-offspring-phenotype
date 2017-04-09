#PBS -S /bin/bash
#PBS -q batch
#PBS -N cuffquant
#PBS -l nodes=1:ppn=10:rjsnode
#PBS -l walltime=480:00:00
#PBS -l mem=10gb

echo "Starting"
cd $PBS_O_WORKDIR
module load anaconda/2.2.0
module load python/2.7.8
module load boost/1.47.0/gcc447
module load samtools/0.1.19
module load cufflinks/2.2.1

#Cuffquant
echo "Running cuffquant"
time cuffquant -p 10 -b ../../ref/bowtie2/Tguttata_v3.2.4.fa -M ../../ref/misc/Tguttata_v3.2.4_mask.gff -u --library-type fr-firststrand ../../ref/transcriptome/Tguttata_v3.2.4.gff tophat/accepted_hits.bam -o cuffquant  
