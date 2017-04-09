#PBS -S /bin/bash
#PBS -q batch
#PBS -N trim_align
#PBS -l nodes=1:ppn=10:rjsnode
#PBS -l walltime=480:00:00
#PBS -l mem=10gb

echo "Starting"
cd $PBS_O_WORKDIR
module load tophat/2.1.1
module load python/2.7.8
module load java/jdk1.8.0_20
module load fastqc/0.11.4
mkdir raw_QC trimmed_QC

#Unpack fastqs
cd ../fastq/rnaseq
for i in *gz
do
gunzip "$i"
done
cat *.fastq > all.fastq

#QC analysis raw data
echo "First QC analysis"
time fastqc --noextract -t 10 -a ../../../../misc/TruSeq.tsv -o ../../rnaseq/raw_QC all.fastq

#Trim data
echo "Trimming data"
python /usr/local/apps/cutadapt/1.9.dev1/bin/cutadapt --trim-n --max-n 0.1 -m 30 -q 20,20 -e 0.1 -O 3 -a file:../../../../misc/TruSeq.fa -o trimmed.fastq all.fastq
rm all.fastq

#QC analysis trimmed data
echo "Second QC analysis"
time fastqc --noextract -t 10 -a ../../../../misc/TruSeq.tsv -o ../../rnaseq/trimmed_QC trimmed.fastq
cd ../../

#Align data
echo "running tophat"
cd rnaseq
time tophat -F 0.1 --transcriptome-index ../../ref/transcriptome/Tguttata_v3.2.4 --no-coverage-search --b2-sensitive --no-mixed --read-realign-edit-dist 0 -i 70 -M -I 500000 -p 10 --library-type fr-firststrand -o tophat ../../ref/bowtie2/Tguttata_v3.2.4 ../fastq/rnaseq/trimmed.fastq
cd ../

#Repackage everything
echo "Wrapping up"
cd fastq/rnaseq
for i in *fastq
do
gzip "$i"
done

echo "Done"
