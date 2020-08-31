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
module load htseq/0.6.1p1
mkdir fastq_QC

#Unpack fastqs
echo "Unpacking fastq files"
cd ../fastq/rnaseq
if [ -f trimmed2.fastq.gz ]
then
  gunzip trimmed2.fastq.gz
  cd ../../
else
  for i in *gz
  do
    gunzip "$i"
  done
  cat *.fastq > raw.fastq

  #QC analysis raw data
  echo "First QC analysis"
  time fastqc --noextract -t 10 -a ../../../../misc/TruSeq.tsv \
  -o ../../rnaseq/fastq_QC raw.fastq

  #Trim data
  echo "Trimming adapters"
  python /usr/local/apps/cutadapt/1.9.dev1/bin/cutadapt --trim-n --max-n 0.1 \
  -m 30 -q 20,20 -e 0.1 -O 3 -a file:../../../../misc/TruSeq.fa \
  -o trimmed1.fastq raw.fastq
  rm raw.fastq

  #QC analysis adapter trimmed data
  echo "Second QC analysis"
  time fastqc --noextract -t 10 -a ../../../../misc/TruSeq.tsv \
  -o ../../rnaseq/fastq_QC trimmed1.fastq

  echo "Second round of trimming"
  time perl /usr/local/apps/prinseq/0.20.4/prinseq-lite.pl -fastq trimmed1.fastq \
  -out_good trimmed2 -out_bad bad -log -trim_left 13 -trim_right 2 \
  -min_len 30 -trim_ns_right 1
  rm trimmed1.fastq bad.fastq

  #QC analysis final trimmed data
  echo "Third QC analysis"
  time fastqc --noextract -t 10 -a ../../../../misc/TruSeq.tsv \
  -o ../../rnaseq/fastq_QC trimmed2.fastq
  cd ../../
fi

#Align data
echo "Running tophat"
cd rnaseq
time tophat -F 0.1 --transcriptome-index ../../ref/transcriptome/Tguttata_v3.2.4 \
--no-coverage-search --b2-very-sensitive --no-mixed --read-realign-edit-dist 0 \
-i 70 -M -I 500000 -p 10 --library-type fr-firststrand -o tophat \
../../ref/bowtie2/Tguttata_v3.2.4 ../fastq/rnaseq/trimmed2.fastq

#Count reads
echo "Counting reads mapping to genes"
htseq-count -f bam -s reverse -t exon -i gene_id -m union \
tophat/accepted_hits.bam ../../ref/transcriptome/Tguttata_v3.2.4.gff \
> gene_counts.tsv

#Repackage everything
echo "Wrapping up"
cd ../fastq/rnaseq
rm trimmed1.fastq.log
for i in *fastq
do
  gzip "$i"
done
cd ../../

echo "Done"
