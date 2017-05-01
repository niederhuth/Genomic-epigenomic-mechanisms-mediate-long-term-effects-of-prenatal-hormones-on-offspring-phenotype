#PBS -S /bin/bash
#PBS -q batch
#PBS -N salmon
#PBS -l nodes=1:ppn=10:rjsnode
#PBS -l walltime=480:00:00
#PBS -l mem=10gb

echo "Unpacking fastq files"
cd ../fastq/rnaseq
if [ -s trimmed2.fastq.gz ]
then
  gunzip trimmed2.fastq.gz
  cd ../../rnaseq
else
  echo "trimmed2.fastq not found"
  exit
fi 

module load salmon/0.8.2
salmon quant -p 10 --posBias --seqBias --gcBias --fldMean 200 --fieldSD 20 \
--numBiasSamples 10000000 --numGibbsSamples 10 \
-i ../../ref/transcriptome/Tguttata_index -l SR \
-r ../fastq/rnaseq/trimmed2.fastq -o salmon_gibbs

salmon quant -p 10 --posBias --seqBias --gcBias --fldMean 200 --fieldSD 20 \
--numBiasSamples 10000000 --numBootstraps 10 \
-i ../../ref/transcriptome/Tguttata_index -l SR \
-r ../fastq/rnaseq/trimmed2.fastq -o salmon_boots 

