#!/bin/bash -login
#PBS -l walltime=96:00:00
#PBS -l nodes=1:ppn=6
#PBS -l mem=36gb
#PBS -N methylpy

cd $PBS_O_WORKDIR

sample=$(pwd | sed s/^.*\\///)

echo "Running $sample"
cd fastq/methylCseq
if ls *fastq.gz >/dev/null 2>&1
then
	echo "Data present"
else
	echo "Downloading from SRA"
	python ../../script/download_fastq.py "$sample"_methylCseq
fi
if ls *sra >/dev/null 2>&1
then
	echo "Running fastq-dump"
	module load SRAToolkit/2.8.2
	for i in *sra
	do
		fastq-dump --split-3 $i
		rm $i
	done
fi
echo "Unpacking fastqs"
for i in *fastq.gz
do
	gunzip $i
done
cd ../../

mkdir methylCseq
cd methylCseq
echo "Running methylpy"
module load SAMTools/1.5
module load bowtie2/2.3.1
module load Java/1.8.0_31
methylpy paired-end-pipeline \
	--read1-files ../fastq/methylCseq/*_R1_*.fastq \
	--read2-files ../fastq/methylCseq/*_R2_*.fastq \
	--libraries "libA" \
	--sample $sample \
	--forward-ref ../../ref/methylCseq/Tguttata_f \
	--reverse-ref ../../ref/methylCseq/Tguttata_r \
	--ref-fasta ../../ref/methylCseq/Tguttata.fa \
	--num-procs 6 \
	--sort-mem 5 \
	--trim-reads True \
	--path-to-cutadapt "" \
	--adapter-seq-read1 AGATCGGAAGAGCACACGTCTGAAC \
	--adapter-seq-read2 AGATCGGAAGAGCGTCGTGTAGGGA \
	--max-adapter-removal 1 \
	--overlap-length 3 \
	--error-rate 0.1 \
	--min-qual-score 10 \
	--min-read-len 30 \
	--bowtie2 True \
	--path-to-aligner "" \
	--aligner-options "" \
	--merge-by-max-mapq True \
	--path-to-samtools "" \
	--remove-clonal True \
	--keep-clonal-stats True \
	--path-to-picard /mnt/home/niederhu/.local/java/lib \
	--java-options "" \
	--binom-test True \
	--num-upstream-bases 0 \
	--num-downstream-bases 2 \
	--unmethylated-control "ChrL" \
	--min-base-quality 1 \
	--min-cov 3 \
	--sig-cutoff .01 \
	--path-to-output "" \
	--keep-temp-files False \
	--generate-allc-file True \
	--generate-mpileup-file True \
	--remove-chr-prefix False \
	--compress-output True

cd ../fastq/methylCseq
for i in *fastq
do
	gzip $i
done


