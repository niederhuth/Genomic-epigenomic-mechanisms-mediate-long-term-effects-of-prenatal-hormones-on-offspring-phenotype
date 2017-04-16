#PBS -S /bin/bash
#PBS -q htseq-count
#PBS -N jobname
#PBS -l nodes=1:ppn=1:rjsnode
#PBS -l walltime=4:00:00
#PBS -l mem=2gb

cd $PBS_O_WORKDIR

module load python/2.7.8
module load htseq/0.6.1p1

htseq-count -f bam -s no -t exon -i gene_id -m union tophat/accepted_hits.bam ../../ref/transcriptome/Tguttata_v3.2.4.gff > gene_counts2.tsv

