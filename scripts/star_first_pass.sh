#!/bin/bash --login
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100GB
#SBATCH --job-name star_first_pass
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR

#set Variables
samples="FCH1 FCH2 FCH3 FCNT1 FCNT2 FCNT3 FTH1 FTH2 FTH3 FTNT1 FTNT2 FTNT3 MCH1 MCH2 MCH3 MCNT1 MCNT2 MCNT3 MTH1 MTH2 MTH3 MTNT1 MTNT2 MTNT3"
fastq=$(for i in $samples; do awk -v a=$i 'BEGIN {print "../"a"/fastq/rnaseq/trimmed.fastq.gz"}'; done | tr '\n' ',' | sed s/,$//)
index="../ref/star"
output="./"

#Align data
echo "Running star"
STAR --runThreadN 20 --runMode alignReads --genomeDir $index --readFilesIn $fastq --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignIntronMin 20 --alignIntronMax 1000000 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 

#filter out MT junctions
awk '$1!="MT"' SJ.out.tab > SJ.filtered.tab

