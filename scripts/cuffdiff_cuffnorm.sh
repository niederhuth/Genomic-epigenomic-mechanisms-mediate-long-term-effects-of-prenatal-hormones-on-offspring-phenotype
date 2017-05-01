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

#cuffdiff all groups
echo "Running cuffdiff on all groups (pooled)"
time cuffdiff -p 20 -L FCH,MCH,FTH,MTH,FCNT,MCNT,FTNT,MTNT \
--dispersion-method pooled --library-norm-method geometric \
--min-reps-for-js-test 3 -b ../ref/bowtie2/Tguttata_v3.2.4.fa -M \
../ref/misc/Tguttata_v3.2.4_mask.gff -u --library-type fr-firststrand \
-o cuffdiff_all ../ref/transcriptome/Tguttata_v3.2.4.gff ../FCH1/rnaseq/cuffquant/abundances.cxb,../FCH2/rnaseq/cuffquant/abundances.cxb,../FCH3/rnaseq/cuffquant/abundances.cxb ../MCH1/rnaseq/cuffquant/abundances.cxb,../MCH2/rnaseq/cuffquant/abundances.cxb,../MCH3/rnaseq/cuffquant/abundances.cxb ../FTH1/rnaseq/cuffquant/abundances.cxb,../FTH2/rnaseq/cuffquant/abundances.cxb,../FTH3/rnaseq/cuffquant/abundances.cxb ../MTH1/rnaseq/cuffquant/abundances.cxb,../MTH2/rnaseq/cuffquant/abundances.cxb,../MTH3/rnaseq/cuffquant/abundances.cxb ../FCNT1/rnaseq/cuffquant/abundances.cxb,../FCNT2/rnaseq/cuffquant/abundances.cxb,../FCNT3/rnaseq/cuffquant/abundances.cxb ../MCNT1/rnaseq/cuffquant/abundances.cxb,../MCNT2/rnaseq/cuffquant/abundances.cxb,../MCNT3/rnaseq/cuffquant/abundances.cxb ../FTNT1/rnaseq/cuffquant/abundances.cxb,../FTNT2/rnaseq/cuffquant/abundances.cxb,../FTNT3/rnaseq/cuffquant/abundances.cxb ../MTNT1/rnaseq/cuffquant/abundances.cxb,../MTNT2/rnaseq/cuffquant/abundances.cxb,../MTNT3/rnaseq/cuffquant/abundances.cxb

#cuffdiff tissue & treatment
echo "Running cuffdiff on tissue and treatment (pooled)"
time cuffdiff -p 20 -L CH,TH,CNT,TNT --dispersion-method pooled \
--library-norm-method geometric --min-reps-for-js-test 6 \
-b ../ref/bowtie2/Tguttata_v3.2.4.fa -M ../ref/misc/Tguttata_v3.2.4_mask.gff \
-u --library-type fr-firststrand \
-o cuffdiff_tissue_treatment ../ref/transcriptome/Tguttata_v3.2.4.gff ../FCH1/rnaseq/cuffquant/abundances.cxb,../FCH2/rnaseq/cuffquant/abundances.cxb,../FCH3/rnaseq/cuffquant/abundances.cxb,../MCH1/rnaseq/cuffquant/abundances.cxb,../MCH2/rnaseq/cuffquant/abundances.cxb,../MCH3/rnaseq/cuffquant/abundances.cxb ../FTH1/rnaseq/cuffquant/abundances.cxb,../FTH2/rnaseq/cuffquant/abundances.cxb,../FTH3/rnaseq/cuffquant/abundances.cxb,../MTH1/rnaseq/cuffquant/abundances.cxb,../MTH2/rnaseq/cuffquant/abundances.cxb,../MTH3/rnaseq/cuffquant/abundances.cxb ../FCNT1/rnaseq/cuffquant/abundances.cxb,../FCNT2/rnaseq/cuffquant/abundances.cxb,../FCNT3/rnaseq/cuffquant/abundances.cxb,../MCNT1/rnaseq/cuffquant/abundances.cxb,../MCNT2/rnaseq/cuffquant/abundances.cxb,../MCNT3/rnaseq/cuffquant/abundances.cxb ../FTNT1/rnaseq/cuffquant/abundances.cxb,../FTNT2/rnaseq/cuffquant/abundances.cxb,../FTNT3/rnaseq/cuffquant/abundances.cxb,../MTNT1/rnaseq/cuffquant/abundances.cxb,../MTNT2/rnaseq/cuffquant/abundances.cxb,../MTNT3/rnaseq/cuffquant/abundances.cxb

#cuffdiff Hypothalamus
echo "Running cuffdiff on Hypothalamus samples (pooled)"
time cuffdiff -p 20 -L FCH,MCH,FTH,MTH, \
--dispersion-method pooled --library-norm-method geometric \
--min-reps-for-js-test 3 -b ../ref/bowtie2/Tguttata_v3.2.4.fa \
-M ../ref/misc/Tguttata_v3.2.4_mask.gff -u --library-type fr-firststrand \
-o cuffdiff_Hy ../ref/transcriptome/Tguttata_v3.2.4.gff ../FCH1/rnaseq/cuffquant/abundances.cxb,../FCH2/rnaseq/cuffquant/abundances.cxb,../FCH3/rnaseq/cuffquant/abundances.cxb ../MCH1/rnaseq/cuffquant/abundances.cxb,../MCH2/rnaseq/cuffquant/abundances.cxb,../MCH3/rnaseq/cuffquant/abundances.cxb ../FTH1/rnaseq/cuffquant/abundances.cxb,../FTH2/rnaseq/cuffquant/abundances.cxb,../FTH3/rnaseq/cuffquant/abundances.cxb ../MTH1/rnaseq/cuffquant/abundances.cxb,../MTH2/rnaseq/cuffquant/abundances.cxb,../MTH3/rnaseq/cuffquant/abundances.cxb

#cuffdiff Nucleus taeniae
echo "Running cuffdiff on Nucleus taeniae samples (pooled)"
time cuffdiff -p 20 -L FCNT,MCNT,FTNT,MTNT -\
-dispersion-method pooled --library-norm-method geometric \
--min-reps-for-js-test 3 -b ../ref/bowtie2/Tguttata_v3.2.4.fa \
-M ../ref/misc/Tguttata_v3.2.4_mask.gff -u --library-type fr-firststrand \
-o cuffdiff_NT ../ref/transcriptome/Tguttata_v3.2.4.gff ../FCNT1/rnaseq/cuffquant/abundances.cxb,../FCNT2/rnaseq/cuffquant/abundances.cxb,../FCNT3/rnaseq/cuffquant/abundances.cxb ../MCNT1/rnaseq/cuffquant/abundances.cxb,../MCNT2/rnaseq/cuffquant/abundances.cxb,../MCNT3/rnaseq/cuffquant/abundances.cxb ../FTNT1/rnaseq/cuffquant/abundances.cxb,../FTNT2/rnaseq/cuffquant/abundances.cxb,../FTNT3/rnaseq/cuffquant/abundances.cxb ../MTNT1/rnaseq/cuffquant/abundances.cxb,../MTNT2/rnaseq/cuffquant/abundances.cxb,../MTNT3/rnaseq/cuffquant/abundances.cxb

#cuffdiff all groups
echo "Running cuffdiff on all groups (per-condition)"
time cuffdiff -p 20 -L FCH,MCH,FTH,MTH,FCNT,MCNT,FTNT,MTNT \
--dispersion-method per-condition --library-norm-method geometric \
--min-reps-for-js-test 3 -b ../ref/bowtie2/Tguttata_v3.2.4.fa -M \
../ref/misc/Tguttata_v3.2.4_mask.gff -u --library-type fr-firststrand \
-o cuffdiff_all2 ../ref/transcriptome/Tguttata_v3.2.4.gff ../FCH1/rnaseq/cuffquant/abundances.cxb,../FCH2/rnaseq/cuffquant/abundances.cxb,../FCH3/rnaseq/cuffquant/abundances.cxb ../MCH1/rnaseq/cuffquant/abundances.cxb,../MCH2/rnaseq/cuffquant/abundances.cxb,../MCH3/rnaseq/cuffquant/abundances.cxb ../FTH1/rnaseq/cuffquant/abundances.cxb,../FTH2/rnaseq/cuffquant/abundances.cxb,../FTH3/rnaseq/cuffquant/abundances.cxb ../MTH1/rnaseq/cuffquant/abundances.cxb,../MTH2/rnaseq/cuffquant/abundances.cxb,../MTH3/rnaseq/cuffquant/abundances.cxb ../FCNT1/rnaseq/cuffquant/abundances.cxb,../FCNT2/rnaseq/cuffquant/abundances.cxb,../FCNT3/rnaseq/cuffquant/abundances.cxb ../MCNT1/rnaseq/cuffquant/abundances.cxb,../MCNT2/rnaseq/cuffquant/abundances.cxb,../MCNT3/rnaseq/cuffquant/abundances.cxb ../FTNT1/rnaseq/cuffquant/abundances.cxb,../FTNT2/rnaseq/cuffquant/abundances.cxb,../FTNT3/rnaseq/cuffquant/abundances.cxb ../MTNT1/rnaseq/cuffquant/abundances.cxb,../MTNT2/rnaseq/cuffquant/abundances.cxb,../MTNT3/rnaseq/cuffquant/abundances.cxb

#cuffdiff tissue & treatment
echo "Running cuffdiff on tissue and treatment (per-condition)"
time cuffdiff -p 20 -L CH,TH,CNT,TNT --dispersion-method per-condition \
--library-norm-method geometric --min-reps-for-js-test 6 \
-b ../ref/bowtie2/Tguttata_v3.2.4.fa -M ../ref/misc/Tguttata_v3.2.4_mask.gff \
-u --library-type fr-firststrand \
-o cuffdiff_tissue_treatment2 ../ref/transcriptome/Tguttata_v3.2.4.gff ../FCH1/rnaseq/cuffquant/abundances.cxb,../FCH2/rnaseq/cuffquant/abundances.cxb,../FCH3/rnaseq/cuffquant/abundances.cxb,../MCH1/rnaseq/cuffquant/abundances.cxb,../MCH2/rnaseq/cuffquant/abundances.cxb,../MCH3/rnaseq/cuffquant/abundances.cxb ../FTH1/rnaseq/cuffquant/abundances.cxb,../FTH2/rnaseq/cuffquant/abundances.cxb,../FTH3/rnaseq/cuffquant/abundances.cxb,../MTH1/rnaseq/cuffquant/abundances.cxb,../MTH2/rnaseq/cuffquant/abundances.cxb,../MTH3/rnaseq/cuffquant/abundances.cxb ../FCNT1/rnaseq/cuffquant/abundances.cxb,../FCNT2/rnaseq/cuffquant/abundances.cxb,../FCNT3/rnaseq/cuffquant/abundances.cxb,../MCNT1/rnaseq/cuffquant/abundances.cxb,../MCNT2/rnaseq/cuffquant/abundances.cxb,../MCNT3/rnaseq/cuffquant/abundances.cxb ../FTNT1/rnaseq/cuffquant/abundances.cxb,../FTNT2/rnaseq/cuffquant/abundances.cxb,../FTNT3/rnaseq/cuffquant/abundances.cxb,../MTNT1/rnaseq/cuffquant/abundances.cxb,../MTNT2/rnaseq/cuffquant/abundances.cxb,../MTNT3/rnaseq/cuffquant/abundances.cxb

#cuffdiff Hypothalamus
echo "Running cuffdiff on Hypothalamus samples (per-condition)"
time cuffdiff -p 20 -L FCH,MCH,FTH,MTH, \
--dispersion-method per-condition --library-norm-method geometric \
--min-reps-for-js-test 3 -b ../ref/bowtie2/Tguttata_v3.2.4.fa \
-M ../ref/misc/Tguttata_v3.2.4_mask.gff -u --library-type fr-firststrand \
-o cuffdiff_Hy2 ../ref/transcriptome/Tguttata_v3.2.4.gff ../FCH1/rnaseq/cuffquant/abundances.cxb,../FCH2/rnaseq/cuffquant/abundances.cxb,../FCH3/rnaseq/cuffquant/abundances.cxb ../MCH1/rnaseq/cuffquant/abundances.cxb,../MCH2/rnaseq/cuffquant/abundances.cxb,../MCH3/rnaseq/cuffquant/abundances.cxb ../FTH1/rnaseq/cuffquant/abundances.cxb,../FTH2/rnaseq/cuffquant/abundances.cxb,../FTH3/rnaseq/cuffquant/abundances.cxb ../MTH1/rnaseq/cuffquant/abundances.cxb,../MTH2/rnaseq/cuffquant/abundances.cxb,../MTH3/rnaseq/cuffquant/abundances.cxb

#cuffdiff Nucleus taeniae
echo "Running cuffdiff on Nucleus taeniae samples (per-condition)"
time cuffdiff -p 20 -L FCNT,MCNT,FTNT,MTNT -\
-dispersion-method per-condition --library-norm-method geometric \
--min-reps-for-js-test 3 -b ../ref/bowtie2/Tguttata_v3.2.4.fa \
-M ../ref/misc/Tguttata_v3.2.4_mask.gff -u --library-type fr-firststrand \
-o cuffdiff_NT2 ../ref/transcriptome/Tguttata_v3.2.4.gff ../FCNT1/rnaseq/cuffquant/abundances.cxb,../FCNT2/rnaseq/cuffquant/abundances.cxb,../FCNT3/rnaseq/cuffquant/abundances.cxb ../MCNT1/rnaseq/cuffquant/abundances.cxb,../MCNT2/rnaseq/cuffquant/abundances.cxb,../MCNT3/rnaseq/cuffquant/abundances.cxb ../FTNT1/rnaseq/cuffquant/abundances.cxb,../FTNT2/rnaseq/cuffquant/abundances.cxb,../FTNT3/rnaseq/cuffquant/abundances.cxb ../MTNT1/rnaseq/cuffquant/abundances.cxb,../MTNT2/rnaseq/cuffquant/abundances.cxb,../MTNT3/rnaseq/cuffquant/abundances.cxb

#cuffnorm all groups
echo "Running cuffnorm on all groups"
time cuffnorm -p 20 -L FCH,MCH,FTH,MTH,FCNT,MCNT,FTNT,MTNT \
--library-norm-method geometric --library-type fr-firststrand -o cuffnorm_all \
../ref/transcriptome/Tguttata_v3.2.4.gff ../FCH1/rnaseq/cuffquant/abundances.cxb,../FCH2/rnaseq/cuffquant/abundances.cxb,../FCH3/rnaseq/cuffquant/abundances.cxb ../MCH1/rnaseq/cuffquant/abundances.cxb,../MCH2/rnaseq/cuffquant/abundances.cxb,../MCH3/rnaseq/cuffquant/abundances.cxb ../FTH1/rnaseq/cuffquant/abundances.cxb,../FTH2/rnaseq/cuffquant/abundances.cxb,../FTH3/rnaseq/cuffquant/abundances.cxb ../MTH1/rnaseq/cuffquant/abundances.cxb,../MTH2/rnaseq/cuffquant/abundances.cxb,../MTH3/rnaseq/cuffquant/abundances.cxb ../FCNT1/rnaseq/cuffquant/abundances.cxb,../FCNT2/rnaseq/cuffquant/abundances.cxb,../FCNT3/rnaseq/cuffquant/abundances.cxb ../MCNT1/rnaseq/cuffquant/abundances.cxb,../MCNT2/rnaseq/cuffquant/abundances.cxb,../MCNT3/rnaseq/cuffquant/abundances.cxb ../FTNT1/rnaseq/cuffquant/abundances.cxb,../FTNT2/rnaseq/cuffquant/abundances.cxb,../FTNT3/rnaseq/cuffquant/abundances.cxb ../MTNT1/rnaseq/cuffquant/abundances.cxb,../MTNT2/rnaseq/cuffquant/abundances.cxb,../MTNT3/rnaseq/cuffquant/abundances.cxb

#cuffnorm tissue & treatment
echo "Running cuffnorm on tissue and treatment"
time cuffnorm -p 20 -L CH,TH,CNT,TNT \
--library-norm-method geometric --library-type fr-firststrand \
-o cuffnorm_tissue_treatment ../ref/transcriptome/Tguttata_v3.2.4.gff ../FCH1/rnaseq/cuffquant/abundances.cxb,../FCH2/rnaseq/cuffquant/abundances.cxb,../FCH3/rnaseq/cuffquant/abundances.cxb,../MCH1/rnaseq/cuffquant/abundances.cxb,../MCH2/rnaseq/cuffquant/abundances.cxb,../MCH3/rnaseq/cuffquant/abundances.cxb ../FTH1/rnaseq/cuffquant/abundances.cxb,../FTH2/rnaseq/cuffquant/abundances.cxb,../FTH3/rnaseq/cuffquant/abundances.cxb,../MTH1/rnaseq/cuffquant/abundances.cxb,../MTH2/rnaseq/cuffquant/abundances.cxb,../MTH3/rnaseq/cuffquant/abundances.cxb ../FCNT1/rnaseq/cuffquant/abundances.cxb,../FCNT2/rnaseq/cuffquant/abundances.cxb,../FCNT3/rnaseq/cuffquant/abundances.cxb,../MCNT1/rnaseq/cuffquant/abundances.cxb,../MCNT2/rnaseq/cuffquant/abundances.cxb,../MCNT3/rnaseq/cuffquant/abundances.cxb ../FTNT1/rnaseq/cuffquant/abundances.cxb,../FTNT2/rnaseq/cuffquant/abundances.cxb,../FTNT3/rnaseq/cuffquant/abundances.cxb,../MTNT1/rnaseq/cuffquant/abundances.cxb,../MTNT2/rnaseq/cuffquant/abundances.cxb,../MTNT3/rnaseq/cuffquant/abundances.cxb

#cuffdiff Hypothalamus
echo "Running cuffnorm on Hypothalamus samples"
time cuffnorm -p 20 -L FCH,MCH,FTH,MTH, \
--library-norm-method geometric--library-type fr-firststrand \
-o cuffnorm_Hy ../ref/transcriptome/Tguttata_v3.2.4.gff ../FCH1/rnaseq/cuffquant/abundances.cxb,../FCH2/rnaseq/cuffquant/abundances.cxb,../FCH3/rnaseq/cuffquant/abundances.cxb ../MCH1/rnaseq/cuffquant/abundances.cxb,../MCH2/rnaseq/cuffquant/abundances.cxb,../MCH3/rnaseq/cuffquant/abundances.cxb ../FTH1/rnaseq/cuffquant/abundances.cxb,../FTH2/rnaseq/cuffquant/abundances.cxb,../FTH3/rnaseq/cuffquant/abundances.cxb ../MTH1/rnaseq/cuffquant/abundances.cxb,../MTH2/rnaseq/cuffquant/abundances.cxb,../MTH3/rnaseq/cuffquant/abundances.cxb

#cuffdiff Nucleus taeniae
echo "Running cuffnorm on Nucleus taeniae samples"
time cuffnorm -p 20 -L FCNT,MCNT,FTNT,MTNT \
--library-norm-method geometric --library-type fr-firststrand -o cuffnorm_NT ../ref/transcriptome/Tguttata_v3.2.4.gff ../FCNT1/rnaseq/cuffquant/abundances.cxb,../FCNT2/rnaseq/cuffquant/abundances.cxb,../FCNT3/rnaseq/cuffquant/abundances.cxb ../MCNT1/rnaseq/cuffquant/abundances.cxb,../MCNT2/rnaseq/cuffquant/abundances.cxb,../MCNT3/rnaseq/cuffquant/abundances.cxb ../FTNT1/rnaseq/cuffquant/abundances.cxb,../FTNT2/rnaseq/cuffquant/abundances.cxb,../FTNT3/rnaseq/cuffquant/abundances.cxb ../MTNT1/rnaseq/cuffquant/abundances.cxb,../MTNT2/rnaseq/cuffquant/abundances.cxb,../MTNT3/rnaseq/cuffquant/abundances.cxb

echo "Done"
