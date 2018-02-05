#!/bin/bash -login
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=6gb
#PBS -N setup

cd $PBS_O_WORKDIR
module load bowtie2/2.3.1
module load SAMTools/1.5

mkdir ref
mkdir ref/sequences ref/methylCseq ref/annotations
cd ref/sequences

#urls
name='Tguttata'
genome='Taeniopygia_guttata.taeGut3.2.4.dna.toplevel.fa'
genome_url='ftp://ftp.ensembl.org/pub/release-91/fasta/taeniopygia_guttata/dna/Taeniopygia_guttata.taeGut3.2.4.dna.toplevel.fa.gz'
gff='Taeniopygia_guttata.taeGut3.2.4.91.gff3'
gff_url='ftp://ftp.ensembl.org/pub/release-91/gff3/taeniopygia_guttata/Taeniopygia_guttata.taeGut3.2.4.91.gff3.gz'

#Download genome
curl $genome_url > $genome.gz
curl $gff_url > $gff
mv $gff ../annotations/

#Prep genome
gunzip $genome.gz -c > $name.fa
samtools faidx $name.fa
gunzip ../../../misc/ChrL.fa.gz -c > ChrL.fa
cat $name.fa ChrL.fa > ../methylCseq/tmp
rm ChrL.fa
gzip $name.fa

#methylCseq index
cd ../methylCseq
python ../../../scripts/fix_fasta.py -i tmp -o $name.fa
rm tmp
samtools faidx $name.fa
methylpy build-reference --input-files $name.fa \
	--output-prefix $name --bowtie2=True

#Format annotations
cd ../annotations
gunzip $gff -c > $name.gff

#grep biotype=processed_pseudogene tmp > tmp2
#fgrep -v -f tmp2 tmp | awk '$1 != "MT"' | awk '$3 == "gene",/###/' > Tguttata_v3.2.4.gff
#awk '$3 == "chromosome" && $1 == "MT"' tmp > Tguttata_v3.2.4_mask.gff
#fgrep -v -f Tguttata_v3.2.4.gff tmp | grep -v \# | awk '$3 != "biological_region" && $3 != "chromosome"' >> Tguttata_v3.2.4_mask.gff
#rm tmp tmp2
#time /usr/local/apps/cufflinks/2.2.1/bin/gffread Tguttata_v3.2.4.gff -T -o Tguttata_v3.2.4.gtf
#time ../../../scripts/cuffdiff_gtf_attributes --input=Tguttata_v3.2.4.gtf --output=Tguttata_v3.2.4_fixed.gtf
#module load TopHat2/2.1.1
#module load cufflinks/2.2.1



#Setup Sample folders
cd ../../
samples=$(sed '1d' ../misc/samples.csv | cut -d ',' -f6 | tr '\n' ' ')
for i in $samples
do
	mkdir $i
done 


