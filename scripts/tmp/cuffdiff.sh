#PBS -S /bin/bash
#PBS -q batch
#PBS -N cuffdiff_cuffnorm
#PBS -l nodes=1:ppn=22:rjsnode
#PBS -l walltime=480:00:00
#PBS -l mem=10gb

echo "Starting"
cd $PBS_O_WORKDIR
module load anaconda/2.2.0
module load python/2.7.8
module load boost/1.47.0/gcc447
module load samtools/0.1.19
module load cufflinks/2.2.1

#cuffdiff
echo "Running cuffdiff on all pair-wise comparisons"
time cuffdiff -p 20 -L FCH,MCH,FTH,MTH,FCNT,MCNT,FTNT,MTNT \
--dispersion-method per-condition --library-norm-method geometric \
--min-reps-for-js-test 3 -b ../ref/bowtie2/Tguttata_v3.2.4.fa -M \
../ref/misc/Tguttata_v3.2.4_mask.gff -u --library-type fr-firststrand \
-o cuffdiff_all2 ../ref/transcriptome/Tguttata_v3.2.4.gff ../FCH1/rnaseq/cuffquant/abundances.cxb,../FCH2/rnaseq/cuffquant/abundances.cxb,../FCH3/rnaseq/cuffquant/abundances.cxb ../MCH1/rnaseq/cuffquant/abundances.cxb,../MCH2/rnaseq/cuffquant/abundances.cxb,../MCH3/rnaseq/cuffquant/abundances.cxb ../FTH1/rnaseq/cuffquant/abundances.cxb,../FTH2/rnaseq/cuffquant/abundances.cxb,../FTH3/rnaseq/cuffquant/abundances.cxb ../MTH1/rnaseq/cuffquant/abundances.cxb,../MTH2/rnaseq/cuffquant/abundances.cxb,../MTH3/rnaseq/cuffquant/abundances.cxb ../FCNT1/rnaseq/cuffquant/abundances.cxb,../FCNT2/rnaseq/cuffquant/abundances.cxb,../FCNT3/rnaseq/cuffquant/abundances.cxb ../MCNT1/rnaseq/cuffquant/abundances.cxb,../MCNT2/rnaseq/cuffquant/abundances.cxb,../MCNT3/rnaseq/cuffquant/abundances.cxb ../FTNT1/rnaseq/cuffquant/abundances.cxb,../FTNT2/rnaseq/cuffquant/abundances.cxb,../FTNT3/rnaseq/cuffquant/abundances.cxb ../MTNT1/rnaseq/cuffquant/abundances.cxb,../MTNT2/rnaseq/cuffquant/abundances.cxb,../MTNT3/rnaseq/cuffquant/abundances.cxb

echo "Done"
