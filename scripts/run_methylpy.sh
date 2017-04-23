#PBS -S /bin/bash
#PBS -q batch
#PBS -N methylpy
#PBS -l nodes=1:ppn=12:HIGHMEM
#PBS -l walltime=480:00:00
#PBS -l mem=100gb

cd $PBS_O_WORKDIR
module load python/2.7.8
echo "Starting"
mkdir allc reports
sample=$(pwd | sed s/.*data\\/// | sed s/\\/methylCseq//)

#Uncompress fastq files
cd ../fastq/methylCseq
echo "Uncompressing fastq files"
for i in *gz
do
  gunzip "$i"
done

#Reverse complement
echo "Reverse complementing the paired end read"
for i in *_R2_001.fastq
do
  output=$(echo "$i" | sed s/.fastq/_rc.fastq/)
  python /usr/local/apps/cutadapt/1.9.dev1/bin/cutadapt -m 30 \
  -a AGATCGGAAGAGCGTCGTGTAGGGA -o tmp.fastq "$i"
  time /usr/local/apps/fastx/0.0.14/bin/fastx_reverse_complement -i tmp.fastq -o "$output"
  rm tmp.fastq
  gzip "$i"
done

#Run methylpy
cd ../../methylCseq
echo "run methylpy"
module load python/2.7.8
python ../../../scripts/run_methylpy.py "$sample" "../fastq/methylCseq/*.fastq" \
"../../ref/methylCseq/Tguttata_v3.2.4" "10" "9" "AGATCGGAAGAGCACACGTCTGAAC" \
"ChrL" > reports/"$sample"_output.txt

#Format allc files
echo "Formatting allc files"
mv allc_* allc/
cd allc
mkdir tmp
head -1 allc_"$sample"_ChrL.tsv > tmp/header
for i in allc_"$sample"_*
do
  sed '1d' "$i" > tmp/"$i"
done

tar -cjvf "$sample"_allc.tar.bz2 allc_"$sample"_*
rm allc_*
cd tmp
rm allc_"$sample"_ChrL.tsv allc_"$sample"_MT.tsv
cat header allc_* > ../"$sample"_allc_total.tsv
cd ../
rm -R tmp
tar -cjvf "$sample"_allc_total.tar.bz2 "$sample"_allc_total.tsv
cd ../

#Cleanup directory
echo "Cleaning up intermediate files"
rm *mpileup* *.bam *.bam.bai

#Compress fastq files
cd ../fastq/methylCseq
echo "Compressing fastq files"
rm *_rc.fastq
for i in *fastq
do
  gzip "$i"
done
cd ../../methylCseq

echo "done"
