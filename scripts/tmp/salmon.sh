#PBS -S /bin/bash
#PBS -q batch
#PBS -N salmon
#PBS -l nodes=1:ppn=10:rjsnode
#PBS -l walltime=480:00:00
#PBS -l mem=10gb

cd $PBS_O_WORKDIR

cd ../fastq/rnaseq
if [ -s trimmed2.fastq.gz ]
then
  echo "Unpacking fastq"
  gunzip trimmed2.fastq.gz
  cd ../../rnaseq
elif [ -s trimmed2.fastq ]
then
  echo "trimmed2.fastq exists"
  cd ../../rnaseq
else
  echo "trimmed2.fastq not found"
  exit
fi 

module load salmon/0.8.2
salmon quant -p 10 --posBias --seqBias --gcBias --fldMean 200 --fldSD 20 \
--numBiasSamples 10000000 --numGibbsSamples 10 \
-i ../../ref/transcriptome/Tguttata_index -l SR \
-r ../fastq/rnaseq/trimmed2.fastq -o salmon_gibbs

salmon quant -p 10 --posBias --seqBias --gcBias --fldMean 200 --fldSD 20 \
--numBiasSamples 10000000 --numBootstraps 10 \
-i ../../ref/transcriptome/Tguttata_index -l SR \
-r ../fastq/rnaseq/trimmed2.fastq -o salmon_boots 

cd ../fastq/rnaseq
gzip trimmed2.fastq
cd ../../rnaseq

