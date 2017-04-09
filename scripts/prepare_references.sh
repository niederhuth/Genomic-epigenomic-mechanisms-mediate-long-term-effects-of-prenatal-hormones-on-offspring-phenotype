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
module load tophat/2.1.1
module load python/2.7.8

#Unpack needed fasta files
cd ../../misc
gunzip ChrL.fa.gz
cd ../data/ref

#Format genome files
echo "Formatting the genome files"
cd bowtie2
gunzip Taeniopygia_guttata.taeGut3.2.4.dna.toplevel.fa.gz
sed s/_random/random/g Taeniopygia_guttata.taeGut3.2.4.dna.toplevel.fa > Tguttata_v3.2.4.fa
gzip Taeniopygia_guttata.taeGut3.2.4.dna.toplevel.fa

echo "Preparing fai, dict, and genome files"
time samtools faidx Tguttata_v3.2.4.fa
cut -f1,2 Tguttata_v3.2.4.fa.fai > Tguttata_v3.2.4.genome
time java -jar /usr/local/apps/picard/1.87/CreateSequenceDictionary.jar R=Tguttata_v3.2.4.fa O=Tguttata_v3.2.4.dict

echo "Making bowtie2 index"
time bowtie2-build Tguttata_v3.2.4.fa Tguttata_v3.2.4
cd ../

#Format annotations
echo "Formatting annotations"
cd misc
gunzip Taeniopygia_guttata.taeGut3.2.4.88.gff3.gz
sed s/_random/random/g Taeniopygia_guttata.taeGut3.2.4.88.gff3 > tmp
gzip Taeniopygia_guttata.taeGut3.2.4.88.gff3
awk '$3 == "gene",/###/' tmp > Tguttata_v3.2.4.gff
awk '$3 == "chromosome" && $1 == "MT",/###/' tmp > Tguttata_v3.2.4_mask.gff
awk '$3 == "miRNA_gene" || $3 == "pseudogene" || $3 == "rRNA_gene" || $3 == "snoRNA_gene" || $3 == "snRNA_gene" || $3 == "mt_gene",/###/' tmp >> Tguttata_v3.2.4_mask.gff
rm tmp
#time ../../../scripts/gffread Tguttata_v3.2.4.gff -T -o Tguttata_v3.2.4.gtf
#time ../../../scripts/cuffdiff_gtf_attributes --input=Tguttata_v3.2.4.gtf --output=Tguttata_v3.2.4_fixed.gtf
cd ../

#Prepare transcriptome index
echo "Making transcriptome index"
cd transcriptome
time tophat2 -G ../misc/Tguttata_v3.2.4_fixed.gtf --transcriptome-index tmp/Tguttata_v3.2.4 ../bowtie2/Tguttata_v3.2.4
mv tmp/* ./
rm -R tophat_out/ tmp/
cd ../

#Prepare methylCseq index
echo "Making methylCseq indexes"
cd methylCseq
cat ../bowtie2/Tguttata_v3.2.4.fa ../../../misc/ChrL.fa > tmp
python ../../../scripts/fix_fasta.py -i tmp -o Tguttata_v3.2.4.fasta
rm tmp
time samtools faidx Tguttata_v3.2.4.fasta
python ../../../scripts/build_Tguttata_methylCseq_index.py
cd ../../

#Wrap up misc files
echo "Wrapping up"
cd ../misc
gzip ChrL.fa
cd ../data
echo "Done"
