#!/bin/bash --login
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=10GB
#SBATCH --job-name trim_reads
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR

#Fastq files
r="raw.fastq.gz"
t="trimmed.fastq.gz"

#Adapter sequences
a="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"

#cat fastq files together
cd fastq/rnaseq
cat *fastq.gz > $r

#Cutadapt
echo "Running cutadapt"
cutadapt -j 10 --trim-n -m 30 -a $a -o $t $r

#Fastqc
mkdir fastqc
echo "Running fastqc"
fastqc -t 10 -o fastqc/ $t $r

echo "Done"
