#PBS -S /bin/bash
#PBS -q batch
#PBS -N build indexes
#PBS -l nodes=1:ppn=2:rjsnode
#PBS -l walltime=480:00:00
#PBS -l mem=10gb

echo "Starting"
cd $PBS_O_WORKDIR
module load bowtie2/2.2.9
module load samtools/1.2
module load picard/2.4.1
module load tophat/2.0.13
module load python/2.7.8

#Format files
echo "Formatting files"
cd bowtie2
gunzip Taeniopygia_guttata.taeGut3.2.4.dna.toplevel.fa.gz
mv Taeniopygia_guttata.taeGut3.2.4.dna.toplevel.fa Tguttata_v3.2.4.fa
sed -i s/_random/random/g Tguttata_v3.2.4.fa
cat Tguttata_v3.2.4.fa ../misc/chrL.fa > ../methylCseq/Tguttata_v3.2.4_all_chr.fa
cd ../misc
gunzip Taeniopygia_guttata.taeGut3.2.4.86.gff3.gz
mv Taeniopygia_guttata.taeGut3.2.4.86.gff3 Tguttata_v3.2.4.gff
sed -i s/_random/random/g Tguttata_v3.2.4.gff
time ../../../scripts/gffread Tguttata_v3.2.4.gff -T -o Tguttata_v3.2.4.gtf
time bash ../../../scripts/cuffdiff_gtf_attributes --input=Tguttata_v3.2.4.gtf --output=Tguttata_v3.2.4_fixed.gtf
sed -i s/_random/random/g Tguttata_v3.2.4_mask.gff
cd ../

#Prepare fai, dict, and genome files
echo "Preparing fai, dict, and genome files"
cd bowtie2
time samtools faidx Tguttata_v3.2.4.fa
time java -jar /usr/local/apps/picard/1.87/CreateSequenceDictionary.jar R=Tguttata_v3.2.4.fa O=Tguttata_v3.2.4.dict
cut -f1,2 Tguttata_v3.2.4.fa.fai > Tguttata_v3.2.4.genome
cd ../methylCseq
time samtools faidx Tguttata_v3.2.4_all_chr.fa
cd ../

#Prepare bowtie2 index
echo "Making bowtie2 index"
cd bowtie2/
time bowtie2-build Tguttata_v3.2.4.fa Tguttata_v3.2.4
cd ../

#Prepare transcriptome index
echo "Making transcriptome index"
cd transcriptome
time tophat2 -G ../misc/Tguttata_v3.2.4_fixed.gtf --transcriptome-index tmp/Tguttata_v3.2.4 ../bowtie2/Tguttata_v3.2.4
mv tmp/* ../transcriptome/
rm -R tophat_out/ tmp/
cd ../

#Prepare methylCseq index
echo "Making methylCseq indexes"
cd methylCseq
python ../../../scripts/build_Tguttata_methylCseq_index.py

echo "Done"
