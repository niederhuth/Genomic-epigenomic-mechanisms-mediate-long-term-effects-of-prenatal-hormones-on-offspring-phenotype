#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20GB
#SBATCH --job-name setup
#SBATCH --output=%x-%j.SLURMout

cd $PBS_O_WORKDIR
PATH="/mnt/home/niederhu/miniconda3/envs/zebra-finch-brains/bin:$PATH"
LD_LIBRARY_PATH="/mnt/home/niederhu/miniconda3/envs/zebra-finch-brains/lib:$LD_LIBRARY_PATH"

#urls
name='Tguttata'
genome='Taeniopygia_guttata.taeGut3.2.4.dna.toplevel.fa'
genome_url='ftp://ftp.ensembl.org/pub/release-91/fasta/taeniopygia_guttata/dna/Taeniopygia_guttata.taeGut3.2.4.dna.toplevel.fa.gz'
gff='Taeniopygia_guttata.taeGut3.2.4.91.gff3.gz'
gff_url='ftp://ftp.ensembl.org/pub/release-91/gff3/taeniopygia_guttata/Taeniopygia_guttata.taeGut3.2.4.91.gff3.gz'

#Download genome
mkdir ref
mkdir ref/sequences ref/methylCseq ref/annotations
cd ref/sequences
curl $genome_url > $genome.gz

#Prep genome
gunzip $genome.gz -c > $name.fa
samtools faidx $name.fa
gunzip ../../../misc/ChrL.fa.gz -c > ChrL.fa
cat $name.fa ChrL.fa > ../methylCseq/tmp
rm ChrL.fa

#methylCseq index
cd ../methylCseq
python ../../../scripts/fix_fasta.py -i tmp -o $name.fa
rm tmp
samtools faidx $name.fa
methylpy build-reference \
	--input-files $name.fa \
	--output-prefix $name \
	--aligner bowtie2

#Format annotations
cd ../annotations
curl $gff_url > $gff
gunzip $gff -c > $name.gff
gffread $name.gff -C -T -o $name.gtf
cd ../

#STAR index
mkdir star
STAR \
	--runThreadN 10 \
	--runMode genomeGenerate \
	--genomeDir star \
	--genomeFastaFiles sequences/"$name".fa \
	--sjdbGTFfile annotations/"$name".gtf \
	--sjdbOverhang 74



#Setup Sample folders
cd ../
samples=$(sed '1d' ../misc/samples.csv | cut -d ',' -f6 | tr '\n' ' ')
for i in $samples
do
	mkdir $i $i/fastq $i/fastq/rnaseq
done 


