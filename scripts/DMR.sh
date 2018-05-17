#!/bin/bash -login
#PBS -l walltime=3:59:00 
#PBS -l nodes=1:ppn=6
#PBS -l mem=10gb
#PBS -N Find_DMRs

#cd $PBS_O_WORKDIR

#SAMPLES="FCH1 FCH2 FCH3 FCNT1 FCNT2 FCNT3 FTH1 FTH2 FTH3 FTNT1 FTNT2 FTNT3 MCH1 MCH2 MCH3 MCNT1 MCNT2 MCNT3 MTH1 MTH2 MTH3 MTNT1 MTNT2 MTNT3"
SAMPLES="FCH1 FTH1"
SAMPLE_CATEGORY=$(seq 1 $(echo $SAMPLES | wc -w) | tr '\n' ',' | sed s/,$//i | sed s/,/,\ /g)
#CHROMS="1 10 11 12 13 14 15 17 18 19 1A 1B 2 20 21 22 23 24 25 26 27 28 3 4 4A 5 6 7 8 9 LG2 LG5 LGE22 Z Un 4random 8random Zrandom 13random 5random 6random 21random 2random 26random 3random 1random 22random 1Arandom 7random 10random 23random 18random 25random LGE22random 9random 15random 2random 20random 11random 4Arandom 14random 17random 27random 19random 28random 16random 24random 1Brandom"
CHROMS="1"
ALLC_FILES=$(for i in $SAMPLES; do awk -v a=$i 'BEGIN {print "../"a"/methylCseq/allc_"a".tsv.gz"}'; done | tr '\n' ' ')
echo $ALLC_FILES

for i in CGN #CHN CNN
do
	echo "Calling $i DMRs"
	methylpy DMRfind --allc-files $ALLC_FILES \
        	         --samples $SAMPLES \
			 --output-prefix $i \
			 --chroms $CHROMS \
                	 --mc-type $i \
                	 --num-procs 6 \
			 --min-cov 1 \
                	 --dmr-max-dist 250 \
                	 --sig-cutoff 0.01 \
			 --num-sims 3000 \
                	 --min-tests 100 \
			 --min-num-dms 0 \
                	 --mc-max-dist 0 \
                	 --keep-temp-files False \
                	 --min-cluster 2 \
			 --seed -1 \
			 --sample-category $SAMPLE_CATEGORY 
			 #--resid-cutoff 0.01 #This is default value, do not change. This is only here to show what the value used in analyses is   
done
