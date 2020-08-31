#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=50GB
#SBATCH --job-name star_second_pass
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR

#set Variables
fastq="../fastq/rnaseq/trimmed.fastq.gz"
index="../../ref/star"
output="./"
junctions="../../rnaseq/SJ.filtered.tab"

#Align data
cd rnaseq
echo "Running star"
STAR --runThreadN 20 --runMode alignReads --genomeDir $index --readFilesIn $fastq --sjdbFileChrStartEnd $junctions --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignIntronMin 20 --alignIntronMax 1000000 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04  --quantMode GeneCounts 

